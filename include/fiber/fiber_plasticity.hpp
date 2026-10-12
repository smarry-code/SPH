#pragma once
// fiber/fiber_plasticity.hpp
// Определяющий закон одноосной связи нити: упруго-пластика, скоростная
// зависимость предела текучести и вязкоупругость (fiber/fiber-plasticity-design.md,
// разделы 1, 3, 6; этап 1 из раздела 10).
//
// ЭТАП 1: чистые функции без зависимости от решателя, FiberState и
// FiberMaterial. Заголовок самодостаточен (только <cmath>/<vector>), поэтому
// тестируется отдельно (tests/test_fiber_plasticity.cpp). На этапе 2
// структуры PlasticParams/RateParams/ViscoParams/FailureParams переедут в
// конец fiber::FiberMaterial, а BondLaw собирается из материала + E.
//
// Реология связи (одноосная, ε = l/l0 - 1):
//   σ     = σ_ep + q + η_v·ε̇
//   σ_ep  = E (ε - ε_p)                      пружина E в серии с пластикой
//   q     : q̇ = E1 ε̇ - q/τ                  ветвь Максвелла
//   η_v·ε̇                                    демпфер Кельвина (в RHS)
// Только растяжение: при ε <= ε_p сила нулевая, σ_total клипуется снизу нулём.
//
// ЕДИНИЦЫ ЭНЕРГИЙ. w_p, w_v -- энергия НА ЕДИНИЦУ ОБЪЁМА связи [Дж/м^3]
// (работа/объём), а не Дж на связь, как в черновике дизайна (раздел 2). Так
// функция не должна знать A и l0; энергия связи = w·A·l0, умножение делает
// вызывающий (Postprocessing этапа 2). Это единственное отклонение
// сигнатуры от дизайна.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace fiber
{

// ----------------------------------------------------------------------
// Параметры (значения по умолчанию = «выключено»)
// ----------------------------------------------------------------------

// Базовый тип 64 бита: иначе -Wpadded (проект собирается с -Weverything)
// ругается на 4 байта выравнивания перед double в PlasticParams/RateParams.
enum class HardeningType : std::int64_t
{
    Linear, // σ_s = sy + H·ε_p        (H = 0 -- идеальная пластичность)
    Power,  // σ_s = sy + K·ε_p^n
    Table   // кусочно-линейная по (table_eps, table_sigma); за последней точкой -- константа
};

struct PlasticParams
{
    HardeningType type = HardeningType::Linear;
    double sy = 0.0; // предел текучести σ_y0 [Па]; sy <= 0 (и пустая таблица) => пластики нет
    double H = 0.0;  // модуль линейного упрочнения [Па]
    double K = 0.0;  // коэффициент степенного упрочнения [Па]
    double n = 1.0;  // показатель степенного упрочнения, 0 < n <= 1
    // Таблица: возрастающие table_eps (первая точка ε_p = 0), неубывающие table_sigma.
    std::vector<double> table_eps;
    std::vector<double> table_sigma;

    bool enabled() const
    {
        if (type == HardeningType::Table) return table_eps.size() >= 2 && table_eps.size() == table_sigma.size();
        return sy > 0.0;
    }
};

enum class RateModel : std::int64_t
{
    None,
    CowperSymonds, // R = 1 + (ε̇_p/C)^{1/p};              a = C [1/с], b = p
    JohnsonCook,   // R = 1 + c·ln(max(ε̇_p/ε̇_0, 1));     a = c,       b = ε̇_0 [1/с]
    Perzyna        // σ_y = σ_s + η_p·ε̇_p^m (аддитивно);  a = η_p [Па·с^m], b = m
};

struct RateParams
{
    RateModel model = RateModel::None;
    double a = 0.0;
    double b = 0.0;
};

struct ViscoParams
{
    double eta = 0.0; // вязкость Кельвина [Па·с]; используется bond_stress
    double xi = 0.0;  // альтернатива: доля критического демпфирования (см. kelvin_eta_from_xi)
    double E1 = 0.0;  // модуль ветви Максвелла [Па]; 0 => ветви нет
    double tau = 0.0; // время релаксации [с]; нужно tau > 0 при E1 > 0
};

struct FailureParams
{
    double eps_p = 0.0;    // предельная ε_p; 0 => критерий выключен
    double sigma_u = 0.0;  // предельное σ_ep; 0 => критерий выключен
};

inline const PlasticParams kNoPlastic{};
inline const RateParams kNoRate{};
inline const ViscoParams kNoVisco{};
inline const FailureParams kNoFailure{};

// Закон связи = E + все блоки. НЕВЛАДЕЮЩЕЕ представление: хранит только ссылки
// (никаких копий таблицы упрочнения, никаких аллокаций), поэтому его можно
// собирать на каждую связь в горячем цикле. Ссылаемые параметры должны жить
// дольше BondLaw (на этапе 2 это поля FiberMaterial).
struct BondLaw
{
    double E; // модуль Юнга [Па]
    const PlasticParams& plastic;
    const RateParams& rate;
    const ViscoParams& visco;
    const FailureParams& failure;

    BondLaw(double E_, const PlasticParams& p = kNoPlastic, const RateParams& r = kNoRate,
            const ViscoParams& v = kNoVisco, const FailureParams& f = kNoFailure)
        : E(E_), plastic(p), rate(r), visco(v), failure(f)
    {
    }
};

// Состояние связи (раздел 2 дизайна).
struct BondState
{
    double eps_p = 0.0;    // пластическая деформация (= накопленная, ε_p >= 0)
    double q = 0.0;        // напряжение ветви Максвелла [Па]
    double eps_prev = 0.0; // ε предыдущего шага
    double w_p = 0.0;      // накопленная пластическая работа, Дж/м^3
    double w_v = 0.0;      // накопленная вязкая диссипация (Максвелл + Кельвин), Дж/м^3
};

// ----------------------------------------------------------------------
// Упрочнение и скоростной множитель
// ----------------------------------------------------------------------

namespace detail
{

// σ_s(ε_p) и dσ_s/dε_p. При n<1 производная в ε_p=0 бесконечна (возвращается inf).
inline void static_yield(const PlasticParams& p, double eps_p, double& s, double& ds)
{
    if (p.type == HardeningType::Linear)
    {
        s = p.sy + p.H * eps_p;
        ds = p.H;
        return;
    }
    if (p.type == HardeningType::Power)
    {
        if (eps_p <= 0.0)
        {
            s = p.sy;
            ds = (p.n < 1.0 && p.K > 0.0) ? std::numeric_limits<double>::infinity() : (p.n < 1.0 ? 0.0 : p.K);
            return;
        }
        s = p.sy + p.K * std::pow(eps_p, p.n);
        ds = p.K * p.n * std::pow(eps_p, p.n - 1.0);
        return;
    }
    // Table
    if (eps_p <= p.table_eps.front())
    {
        s = p.table_sigma.front();
        ds = (p.table_sigma[1] - p.table_sigma[0]) / (p.table_eps[1] - p.table_eps[0]);
        return;
    }
    if (eps_p >= p.table_eps.back())
    {
        s = p.table_sigma.back();
        ds = 0.0;
        return;
    }
    const auto it = std::upper_bound(p.table_eps.begin(), p.table_eps.end(), eps_p);
    const std::size_t k = static_cast<std::size_t>(it - p.table_eps.begin()); // eps[k-1] <= eps_p < eps[k]
    const double slope = (p.table_sigma[k] - p.table_sigma[k - 1]) / (p.table_eps[k] - p.table_eps[k - 1]);
    s = p.table_sigma[k - 1] + slope * (eps_p - p.table_eps[k - 1]);
    ds = slope;
}

} // namespace detail

// Статический предел текучести σ_s(ε_p).
inline double static_yield_stress(const PlasticParams& p, double eps_p)
{
    double s = 0.0, ds = 0.0;
    detail::static_yield(p, eps_p, s, ds);
    return s;
}

// Динамический предел текучести σ_y(ε_p, ε̇_p) = σ_s(ε_p)·R(ε̇_p) (Perzyna: σ_s + η ε̇^m).
// d_deps_p -- (необязательно) частная производная по ε_p, d_dedot -- по ε̇_p
// (inf допустим: тогда Ньютон в решателе уходит в деление отрезка).
inline double plastic_yield_stress(const PlasticParams& p, const RateParams& r, double eps_p, double edot_p,
                                   double* d_deps_p = nullptr, double* d_dedot = nullptr)
{
    double s = 0.0, ds = 0.0;
    detail::static_yield(p, eps_p, s, ds);
    const double inf = std::numeric_limits<double>::infinity();
    const double ed = std::max(edot_p, 0.0);
    double sy = s, dsde = ds, dsdr = 0.0;
    if (r.model == RateModel::CowperSymonds)
    {
        const double ex = 1.0 / r.b;
        const double ratio = ed / r.a;
        const double pw = std::pow(ratio, ex);
        sy = s * (1.0 + pw);
        dsde = ds * (1.0 + pw);
        if (ed > 0.0) dsdr = s * ex * pw / ed;
        else if (ex < 1.0) dsdr = inf;
        else if (ex > 1.0) dsdr = 0.0;
        else dsdr = s / r.a;
    }
    else if (r.model == RateModel::JohnsonCook)
    {
        const double x = ed / r.b;
        if (x > 1.0)
        {
            const double R = 1.0 + r.a * std::log(x);
            sy = s * R;
            dsde = ds * R;
            dsdr = s * r.a / ed;
        }
        else
        {
            sy = s;
            dsde = ds;
            dsdr = 0.0;
        }
    }
    else if (r.model == RateModel::Perzyna)
    {
        const double pw = std::pow(ed, r.b);
        sy = s + r.a * pw;
        dsde = ds;
        if (ed > 0.0) dsdr = r.a * r.b * pw / ed;
        else if (r.b < 1.0) dsdr = inf;
        else if (r.b > 1.0) dsdr = 0.0;
        else dsdr = r.a;
    }
    if (d_deps_p) *d_deps_p = dsde;
    if (d_dedot) *d_dedot = dsdr;
    return sy;
}

// Кельвин через долю критического демпфирования (раздел 1.3):
//   c = η_v·A/l0 = 2ξ·sqrt(k·m_red),  k = E_inst·A/l0  =>  η_v = c·l0/A.
inline double kelvin_eta_from_xi(double xi, double E_inst, double A, double l0, double m_red)
{
    const double k = E_inst * A / l0;
    const double c = 2.0 * xi * std::sqrt(k * m_red);
    return c * l0 / A;
}

// ----------------------------------------------------------------------
// Напряжение связи (раздел 3) -- то, что вызовет RHS-ядро. Только читает.
//   ε <= ε_p        => 0 (связь провисла)
//   иначе σ = max(E(ε-ε_p) + q + η_v·ε̇, 0);  сила = A·σ.
// ----------------------------------------------------------------------
inline double bond_stress(const BondLaw& law, const BondState& s, double eps, double eps_dot)
{
    if (!(eps > s.eps_p)) return 0.0;
    // eta == 0: не умножаем (0*inf = NaN), результат побитово равен упругому закону.
    const double visc = law.visco.eta != 0.0 ? law.visco.eta * eps_dot : 0.0;
    const double sigma = law.E * (eps - s.eps_p) + s.q + visc;
    return std::max(sigma, 0.0);
}

// ----------------------------------------------------------------------
// Конститутивное обновление за шаг (разделы 1.2, 1.3). Чистая функция.
//
// Вход: состояние на начало шага (eps_prev -- ε на начало), ε на конец шага,
// dt. Выход: состояние на конец шага (eps_prev = eps_new). При dt <= 0 (или
// нечисловом) состояние возвращается без изменений.
//
// Пластика (если plastic.enabled()): σ_tr = E(ε - ε_p). Если σ_tr <= σ_s(ε_p)
// -- упруго. Иначе корень g(Δ) = σ_tr - EΔ - σ_y(ε_p+Δ, Δ/dt) = 0 на скобке
// [0, σ_tr/E] (g(0)>0, g(σ_tr/E)<0, g строго убывает => корень единственный):
//   * линейное упрочнение без скорости -- аналитически Δ = (σ_tr-σ_s)/(E+H);
//   * иначе Ньютон с защитой: шаг вне скобки (или производная не конечна)
//     заменяется делением отрезка. Остановка по |g| <= 1e-14·σ_tr или
//     схлопыванию скобки (4·DBL_EPSILON относительно).
// Работа w_p += (σ_start + σ_end)/2·Δ (трапеция; точна для линейного
// упрочнения без скорости; дизайн предлагал верхнюю оценку σ_y(конец)·Δ).
// σ_start -- фактическое напряжение в начале шага, E(ε_prev-ε_p) в пределах
// [σ_s(ε_p), σ_end] (при скоростной зависимости оно выше статического),
// σ_end = σ_tr - EΔ (на корне равно σ_y; так работа ограничена сверху σ_tr·Δ
// и не раздувается при экстремально малом dt).
//
// СОЗНАТЕЛЬНОЕ УТОЧНЕНИЕ ДИЗАЙНА (раздел 1.3, провисшая связь): связь только на
// растяжение, при ε <= ε_p все элементы разгружены. Ветви Максвелла и Кельвина
// видят не полное Δε, а приращение деформации натянутой части
//   Δε_eff = max(ε_new, ε_p_new) - max(ε_prev, ε_p_old).
// Если связь провисла в начале и в конце шага, Δε_eff = 0: q только
// релаксирует экспонентой, диссипация Кельвина не копится. Переход через
// границу ε = ε_p непрерывен (учитывается только натянутая часть шага); при
// пластическом течении ε_new >= ε_p_new, и Δε_eff = ε_new - max(ε_prev, ε_p_old).
//
// Максвелл (E1>0, tau>0): точная экспонента при линейном ε(t) на шаге,
//   q' = e^{-dt/τ} q + E1·Δε_eff·(τ/dt)(1 - e^{-dt/τ}).
// Вязкая диссипация: ∫ q²/(E1 τ) dt в замкнутой форме (на шаге, тот же
// линейный ход) + Кельвин η_v·Δε_eff²/dt. Пластика и Максвелл независимы
// (q не зависит от ε_p).
// ----------------------------------------------------------------------
inline BondState bond_constitutive_update(const BondLaw& law, const BondState& in, double eps_new, double dt)
{
    if (!(dt > 0.0)) return in;
    BondState out = in;
    const double E = law.E;

    // --- пластический возврат ---
    if (law.plastic.enabled())
    {
        const PlasticParams& pl = law.plastic;
        const RateParams& rt = law.rate;
        const double sig_tr = E * (eps_new - in.eps_p);
        const double sig_s0 = static_yield_stress(pl, in.eps_p);
        if (sig_tr > sig_s0)
        {
            double delta = 0.0;
            const double hi0 = sig_tr / E;
            const bool analytic = pl.type == HardeningType::Linear && rt.model == RateModel::None;
            if (analytic)
            {
                delta = (sig_tr - sig_s0) / (E + pl.H);
            }
            else
            {
                double lo = 0.0, hi = hi0;
                double x = 0.5 * hi;
                const double tol = 1e-14 * sig_tr;
                for (int it = 0; it < 200; ++it)
                {
                    double dse = 0.0, dsr = 0.0;
                    const double sy = plastic_yield_stress(pl, rt, in.eps_p + x, x / dt, &dse, &dsr);
                    const double g = sig_tr - E * x - sy;
                    if (std::fabs(g) <= tol) break;
                    if (g > 0.0) lo = x;
                    else hi = x;
                    if (hi - lo <= 4.0 * std::numeric_limits<double>::epsilon() * hi) break;
                    const double dg = -(E + dse + dsr / dt); // g' < 0
                    double xn = x - g / dg;
                    if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
                    x = xn;
                }
                delta = x;
            }
            delta = std::min(std::max(delta, 0.0), hi0);
            const double sig_end = std::min(std::max(sig_tr - E * delta, 0.0), sig_tr);
            const double sig_start = std::min(std::max(E * (in.eps_prev - in.eps_p), sig_s0), sig_end);
            out.eps_p = in.eps_p + delta;
            out.w_p = in.w_p + 0.5 * (sig_start + sig_end) * delta;
        }
    }

    // --- натянутая часть шага (провисшая связь) ---
    const double deps_eff = std::max(eps_new, out.eps_p) - std::max(in.eps_prev, in.eps_p);

    // --- ветвь Максвелла ---
    const ViscoParams& vs = law.visco;
    if (vs.E1 > 0.0 && vs.tau > 0.0)
    {
        const double r = dt / vs.tau;
        const double ex = std::exp(-r);
        const double one_m = -std::expm1(-r); // 1 - e^{-r}, без потери точности при малом r
        // (1 - e^{-r})/r; ряд при малом r (в т.ч. r = 0 при экстремально малом dt: без 0/0)
        const double g1 = r < 1e-8 ? 1.0 - 0.5 * r + r * r / 6.0 : one_m / r;
        out.q = ex * in.q + vs.E1 * deps_eff * g1;
        // q(t) = a e^{-t/τ} + b (1 - e^{-t/τ}), a = q_in, b = E1·ε̇_eff·τ = E1·Δε_eff/r.
        // Интеграл q² на шаге (s = e^{-r}):
        //   a²τ(1-s²)/2 + ab·τ(1-s)² + b²·τ·φ(r),  φ(r) = r - 2(1-e^{-r}) + (1-e^{-2r})/2.
        // Через Δε_eff, без деления на dt (иначе inf*0 при dt ~ 1e-300):
        //   ab(1-s)² = a·E1·Δε_eff·(1-s)·g1,   b²φ = E1²Δε_eff²·ψ,  ψ = φ/r².
        // ψ ~ r/3 при малом r: прямое вычисление теряет точность, поэтому ряд.
        double psi;
        if (r < 0.02)
        {
            psi = r * (1.0 / 3.0 - r * (1.0 / 4.0 - r * (7.0 / 60.0 - r * (1.0 / 24.0 - r * (31.0 / 2520.0 - r / 320.0)))));
        }
        else
        {
            psi = (r - 2.0 * one_m - 0.5 * std::expm1(-2.0 * r)) / (r * r);
        }
        const double a = in.q;
        const double E1d = vs.E1 * deps_eff;
        const double int_q2 =
            vs.tau * (a * a * 0.5 * (-std::expm1(-2.0 * r)) + a * E1d * one_m * g1 + E1d * E1d * psi);
        out.w_v += int_q2 / (vs.E1 * vs.tau);
    }
    if (vs.eta > 0.0) out.w_v += vs.eta * deps_eff * (deps_eff / dt);

    out.eps_prev = eps_new;
    return out;
}

} // namespace fiber
