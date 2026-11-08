#pragma once
// energy_log.hpp
// Диагностический журнал энергетического баланса (CSV) для SPH+Fiber решателя.
//
// Назначение: на выбранных шагах считать запасы энергии системы (кинетическая
// Solid/Fiber, внутренняя, упругая энергия связей, потенциал Герца) и
// накопленные каналы потерь/источников, а затем писать строку CSV с невязкой
// баланса. Нужен для проверки того, куда уходит энергия (трение, XSPH,
// Solid-Solid контакт, разрыв связей, зажим по u).
//
// ВЕСЬ файл находится под #ifdef SPH_ENERGY_LOG: без этого макроса заголовок
// пуст и на сборку/производительность не влияет.
//
// Масса Solid-частицы во всех суммах и каналах -- energy_mass (particle_system.hpp).
// В 2D осесимметрии это масса кольца 2*pi*r0*m: «плоская» масса m = rho0*dx^2
// не даёт сохраняющейся энергии (johnson-cook/03-pmin-energy.md, §2.1), а
// 2*pi*r0*m -- физическая масса, и столбцы Solid -- энергия всего тела (Дж).
// В 1D и 3D energy_mass = m, лог прежний.
//
// Столбцы CSV (энергии в единицах решателя, Дж):
//   step, time, dt   -- номер шага, время, шаг по времени
//   KE_solid         -- sum 1/2 m |v|^2 по реальным частицам Solid
//   KE_fiber         -- sum 1/2 m |v|^2 по бусинам Fiber
//   U_int            -- sum m*u (внутренняя энергия Solid; УЖЕ включает
//                       U_shear_el, Q_av и W_plastic)
//   U_shear_el       -- упругая энергия сдвига sum m*(s:s)/(4 G rho);
//                       ИНФОРМАЦИОННЫЙ, входит в U_int, в E_mech не прибавляется
//   U_bond           -- запасённая энергия связей нитей. Упругий материал:
//                       sum 1/2 (EA/l0)(l-l0)^2 по растянутым связям (l>l0).
//                       Неупругий (PLASTIC/VISCO): 1/2 E (eps-eps_p)^2 A l0 при
//                       eps>eps_p плюс 1/2 q^2/E1 A l0 при E1>0 (пружина
//                       Максвелла); w_p и w_v в запас не входят (см. W_p_fiber)
//   U_hertz_ff       -- потенциал Герца Fiber-Fiber, sum 2/5 F delta (F=k delta^1.5);
//                       при fiber_fiber_contact_type = 0 (Герц) и 4 (Герц + демпфер), иначе 0
//   U_hertz_fs       -- потенциал Герца Fiber-Solid, то же; при
//                       fiber_solid_contact_type = 0 и 4, иначе 0
//   E_mech           -- KE_solid+KE_fiber+U_int+U_bond+U_hertz_ff+U_hertz_fs
//   D_friction       -- накопленная диссипация трения (потери, > 0)
//   D_xsph           -- накопленная диссипация XSPH (потери, > 0); при
//                       xsph_mode = position равна 0 (av не входит в скорость)
//   D_contact_ss     -- накопленная диссипация Solid-Solid контакта (потери, > 0)
//                       (эти три -- работа вклада на средней скорости шага,
//                       см. work_mid_fiber/work_mid_solid)
//   D_break          -- накопленная запасённая энергия порвавшихся связей (> 0; для неупругих --
//                       1/2 E(eps-eps_p)^2 A l0 + пружина Максвелла 1/2 q^2/E1 A l0 по состоянию
//                       после шага разрыва; w_p, w_v этого шага уже в W_p_fiber/W_v_fiber)
//   S_u_clamp        -- накопленный ИСТОЧНИК энергии от зажима u (> 0)
//   Q_av             -- накопленное искусственное тепло (информационно, уже в U_int)
//   W_plastic        -- накопленная работа пластичности (информационно, уже в U_int)
//   residual         -- (E_mech - E_mech0) + D_friction + D_xsph + D_contact_ss
//                       + D_break - S_u_clamp + W_p_fiber + W_v_fiber
//                       + D_break_contact - clip_loss + D_contact_fiber;
//                       E_mech0 берётся с ПЕРВОЙ записанной строки. Для идеально
//                       замкнутого баланса residual ~ 0.
//   W_p_fiber        -- накопленная пластическая работа связей нитей (потери, > 0):
//                       sum по связям и шагам приращений w_p*A*l0, включая связи,
//                       порвавшиеся на этом же шаге (приращение берётся до удаления).
//                       Не путать с W_plastic (Мизес, Solid, справочно).
//   W_v_fiber        -- накопленная вязкая диссипация связей нитей (потери, >= 0):
//                       sum приращений w_v*A*l0 (ветвь Максвелла + демпфер Кельвина).
//                       У упругих нитей оба канала нулевые.
//   D_break_contact  -- знаковая поправка топологии: (U_Герца до - U_Герца после) удаления
//                       порвавшихся связей, накапливается. Сегмент порванной связи пропадает
//                       из контактного поиска (а object_id j-стороны меняется и перестраивает
//                       фильтр соседей), поэтому потенциал Герца Fiber-Fiber/Fiber-Solid
//                       исчезает (или появляется) без работы. На E5 «Максвелл + разрывы»
//                       это главный источник невязки: 0.057 % E0 против D_break 0.014 % E0
//                       (невязка падает с 0.077 % до 0.020 % E0). Может быть < 0.
//   clip_loss        -- поправка отсечки sigma >= 0 (знаковая, обычно > 0 и мала): при
//                       eps > eps_p и sigma_модели = E(eps-eps_p) + q + eta*eps_dot < 0
//                       (быстрая разгрузка Максвелла/Кельвина) сила связи нулевая, но
//                       внутренние ветви (запас, w_v) всё равно меняются на sigma_модели*d(eps)*A*l0;
//                       накапливается эта величина (RHS-ядро сообщает sigma_модели, проход
//                       Postprocessing умножает на d(eps)); в residual вычитается. Не зависит
//                       от dt (E5: ~0.001 % E0, ~1 % от W_v).
//   D_contact_fiber  -- накопленная работа (со знаком минус) нормальной силы контактов нитей,
//                       закон которых НЕ Герц (fiber_*_contact_type = 1..3: pinball, Swegle,
//                       импульсный). Эти силы не потенциальные: запаса энергии у них нет, вся их
//                       работа -- диссипация (> 0: сила действует только на сближении). У закона 4
//                       (Герц + демпфер) упругая часть потенциальна (в U_hertz_*), сюда идёт только
//                       работа демпфера (отдельный проход ядра с damper_only). Работа --
//                       на средней скорости шага (work_mid_fiber/work_mid_solid), по ускорению
//                       контакта после проекции зажимов. При законе Герц канал нулевой, энергия
//                       контакта -- в U_hertz_*. Для разных законов у двух контактов каждый
//                       контакт учитывается своим способом.
//   (W_p_fiber, W_v_fiber, D_break_contact, clip_loss, D_contact_fiber -- последние столбцы:
//    старые CSV без них читаются как нулевые каналы.)
//
// Остаток невязки у нитей с вязкостью/пластикой -- ошибка схемы O(dt) (силы по состоянию на
// начало шага, Кельвин по скорости предиктора, запас связи по x^{n+1}); на E5 с разрывами
// max|residual| сходится ~ пропорционально dt: 2.0e-4 E0 (dt лога), 1.0e-4 (dt/2), 4.3e-5 (dt/8),
// 1.1e-5 (dt/32), 5.3e-6 (dt/64). Подробности и числа: fiber/fiber-plasticity-design.md, раздел 9, риск 4.
//
// Знаки: потери D_*, W_p_fiber, W_v_fiber положительны (энергия уходит из E_mech),
// источник S_u_clamp положителен (энергия добавляется в систему), поэтому он
// вычитается.
//
// Потенциал Герца считается собственным поиском контактов (EnergyStockWorkspace),
// который перестраивается ТОЛЬКО на записываемых строках, а не каждый шаг;
// память поиска переиспользуется между строками.
#ifdef SPH_ENERGY_LOG

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "fiber/fiber_system.hpp"
#include "interactions/fiber_bond.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/fiber_solid_contact.hpp"
#include "materials.hpp"
#include "ode/ode_solver.hpp"
#include "interactions/postprocessing.hpp"
#include "param.hpp"
#include "particle_system.hpp"
#include "rhs/cpu_executor.hpp"
#include "simulation_config.hpp"

namespace energy
{

// Накопленные каналы (потери > 0; u_clamp -- источник > 0;
// q_av, w_plastic -- информационные, уже внутри U_int;
// w_p_fiber, w_v_fiber -- потери нитей (пластика и вязкость связей), входят в баланс).
struct EnergyChannels
{
    double friction = 0, xsph = 0, contact_ss = 0, brk = 0, u_clamp = 0, q_av = 0, w_plastic = 0;
    double w_p_fiber = 0, w_v_fiber = 0;
    double brk_contact = 0, clip = 0; // знаковые поправки учёта (см. шапку)
    double contact_fiber = 0;         // диссипация в контактах нитей с законом 1-3 (см. шапку)

    // Чекпоинт (checkpoint::EnergyState::chan): порядок полей фиксирован -- 12 значений.
    static constexpr std::size_t kCount = 12;
    void to_array(double* a) const
    {
        const double v[kCount] = {friction, xsph, contact_ss, brk, u_clamp, q_av, w_plastic,
                                  w_p_fiber, w_v_fiber, brk_contact, clip, contact_fiber};
        for (std::size_t i = 0; i < kCount; ++i) a[i] = v[i];
    }
    void from_array(const double* a)
    {
        double* f[kCount] = {&friction, &xsph, &contact_ss, &brk, &u_clamp, &q_av, &w_plastic,
                             &w_p_fiber, &w_v_fiber, &brk_contact, &clip, &contact_fiber};
        for (std::size_t i = 0; i < kCount; ++i) *f[i] = a[i];
    }
};

// Мгновенные запасы энергии.
struct EnergyStocks
{
    double ke_solid = 0, ke_fiber = 0, u_int = 0, u_shear_el = 0, u_bond = 0, u_hertz_ff = 0, u_hertz_fs = 0;
    // u_shear_el входит в u_int
    double e_mech() const { return ke_solid + ke_fiber + u_int + u_bond + u_hertz_ff + u_hertz_fs; }
};

// Рабочая память поиска контактов для подсчёта потенциала Герца (живёт в
// логгере, память только растёт).
struct EnergyStockWorkspace
{
    interactions::FiberContactSearchResult ff;
    interactions::FiberSolidContactSearchResult fs;
};

namespace detail
{
struct SolidAcc
{
    double ke = 0.0, u_int = 0.0, u_shear_el = 0.0;
};
} // namespace detail

// Потенциал Герца Fiber-Fiber: U = 2/5 F delta (F = k delta^1.5), и Fiber-Solid, то же.
// Собственный поиск контактов по ТЕКУЩЕЙ топологии (object_id, topology.bonds):
// потенциал зависит от топологии, поэтому при разрыве связи его считают до и после
// (см. FiberEnergyTap::contact_potential и D_break_contact).
inline void compute_hertz_potential(const rhs::CpuExecutor& exec, const SimulationConfig& cfg,
                                    const ParticleSystem& sys, std::size_t ntotal, const fiber::FiberState& fb,
                                    const fiber::FiberTopology& topology,
                                    const std::vector<fiber::FiberMaterial>& materials, EnergyStockWorkspace& ws,
                                    double& u_ff, double& u_fs)
{
    u_ff = 0.0;
    u_fs = 0.0;
    auto sum = [](double a, double b) { return a + b; };
    // Потенциал есть только у закона Герца (тип 0). Законы 1-3 не потенциальные:
    // их работа идёт в канал D_contact_fiber, а не в запас.
    const bool hertz_ff = interactions::fiber_contact_law_has_hertz_potential(cfg.fiber_fiber_contact_type);
    const bool hertz_fs = interactions::fiber_contact_law_has_hertz_potential(cfg.fiber_solid_contact_type);
    if (hertz_ff) interactions::build_fiber_contact_search(exec, fb, topology, materials, ws.ff);
    if (hertz_ff && !ws.ff.primitives.empty())
    {
        u_ff = exec.parallel_reduce(
            ws.ff.primitives.size(), 0.0,
            [&](std::size_t k) -> double {
                const auto geom = interactions::compute_fiber_contact_geometry(fb, materials, ws.ff.primitives[k]);
                if (geom.rejected_by_ownership || geom.degenerate || !geom.contact) return 0.0;
                return 0.4 * geom.F * geom.delta;
            },
            sum);
    }

    if (hertz_fs && ntotal > 0)
    {
        interactions::build_fiber_solid_contact_search(exec, cfg.mats, sys, ntotal, fb, topology, materials, ws.fs);
        if (!ws.fs.primitives.empty())
        {
            u_fs = exec.parallel_reduce(
                ws.fs.primitives.size(), 0.0,
                [&](std::size_t k) -> double {
                    const auto geom = interactions::compute_fiber_solid_contact_geometry(
                        cfg.mats, sys, fb, materials, ws.fs.primitives[k]);
                    if (!geom.contact || !(geom.delta > 0.0)) return 0.0;
                    return 0.4 * geom.F * geom.delta;
                },
                sum);
        }
    }
}

inline EnergyStocks compute_energy_stocks(const rhs::CpuExecutor& exec, const SimulationConfig& cfg,
                                          const ParticleSystem& sys, std::size_t ntotal,
                                          const fiber::FiberState& fb, const fiber::FiberTopology& topology,
                                          const std::vector<fiber::FiberMaterial>& materials,
                                          EnergyStockWorkspace& ws)
{
    EnergyStocks s;

    // --- Solid: KE, внутренняя энергия, упругая энергия сдвига ---
    if (ntotal > 0)
    {
        const detail::SolidAcc acc = exec.parallel_reduce(
            ntotal, detail::SolidAcc{},
            [&](std::size_t i) -> detail::SolidAcc {
                detail::SolidAcc a;
                const double m = energy_mass(sys, i); // 2D осесимметрия: масса кольца
                double v2 = 0.0;
                for (std::size_t d = 0; d < dim; ++d) v2 += sys.vx[d][i] * sys.vx[d][i];
                a.ke = 0.5 * m * v2;
                a.u_int = m * sys.u[i];
                const double G = shear_modulus_by_itype(cfg.mats, sys.itype[i]);
                const double rho = sys.rho[i];
                if (G > 0.0 && rho > 0.0)
                {
                    const double sxx = sys.sxx[i], syy = sys.syy[i], szz = sys.szz[i];
                    const double sxy = sys.sxy[i], syz = sys.syz[i], sxz = sys.sxz[i];
                    a.u_shear_el = m * (sxx * sxx + syy * syy + szz * szz +
                                        2.0 * (sxy * sxy + syz * syz + sxz * sxz)) /
                                   (4.0 * G * rho);
                }
                return a;
            },
            [](detail::SolidAcc a, const detail::SolidAcc& b) {
                a.ke += b.ke;
                a.u_int += b.u_int;
                a.u_shear_el += b.u_shear_el;
                return a;
            });
        s.ke_solid = acc.ke;
        s.u_int = acc.u_int;
        s.u_shear_el = acc.u_shear_el;
    }

    const std::size_t nf = fb.size();
    if (nf == 0) return s;

    auto sum = [](double a, double b) { return a + b; };

    // --- Fiber: кинетическая энергия бусин ---
    s.ke_fiber = exec.parallel_reduce(
        nf, 0.0,
        [&](std::size_t i) -> double {
            double v2 = 0.0;
            for (std::size_t d = 0; d < dim; ++d) v2 += fb.vx[d][i] * fb.vx[d][i];
            return 0.5 * fb.mass[i] * v2;
        },
        sum);

    // --- Fiber: запасённая энергия связей (только растяжение; неупругие -- с вычетом
    //     eps_p и пружиной Максвелла, interactions::fiber_bond_stored_energy) ---
    const std::size_t nb = topology.bonds.size();
    if (nb > 0)
    {
        s.u_bond = exec.parallel_reduce(
            nb, 0.0,
            [&](std::size_t k) -> double {
                const fiber::FiberBond& b = topology.bonds[k];
                const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(fb.material_id[b.i])];
                double l2 = 0.0;
                for (std::size_t d = 0; d < dim; ++d)
                {
                    const double r = fb.x[d][b.j] - fb.x[d][b.i];
                    l2 += r * r;
                }
                return interactions::fiber_bond_stored_energy(mat, b, std::sqrt(l2));
            },
            sum);
    }

    // --- Потенциал Герца Fiber-Fiber и Fiber-Solid ---
    compute_hertz_potential(exec, cfg, sys, ntotal, fb, topology, materials, ws, s.u_hertz_ff, s.u_hertz_fs);
    return s;
}

// Работа одного вклада в скорость за шаг -- на СРЕДНЕЙ скорости шага.
// Обновление решателя: v_new = base + h*a + av (a -- полное ускорение после
// всех стадий, av -- XSPH). Вклад X меняет скорость на dv_X (h*a_X или сам av),
// sum_X dv_X = v_new - base, поэтому
//     W_X = sum_i m_i dv_X,i . v_mid,i,   v_mid = base + (h*a + av)/2,
// и sum_X W_X = KE(v_new) - KE(base) точно. Разбиение не зависит от порядка
// вкладов: бусина, которую трение удерживает против натянутой связи (a_тр = -a_св,
// base = 0), не движется, и обе работы равны нулю. При разбиении «трение после
// остальных сил» трению приписывалась бы ложная потеря 1/2 m h^2 a_св^2 на шаг.

// Fiber: dv_X = h*acc, v_mid = base + h*dv/2 (у нитей нет XSPH).
inline double work_mid_fiber(const rhs::CpuExecutor& exec, const fiber::FiberState& fb,
                             const ode::FiberVelocityUpdate& upd, const std::vector<double> (&dv)[dim],
                             const std::vector<double> (&acc)[dim], std::size_t n)
{
    return exec.parallel_reduce(
        n, 0.0,
        [&](std::size_t i) -> double {
            double w = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double v_mid = (*upd.base[d])[i] + 0.5 * upd.h * dv[d][i];
                w += upd.h * acc[d][i] * v_mid;
            }
            return fb.mass[i] * w;
        },
        [](double a, double b) { return a + b; });
}

// Solid: dv_X = h*inc (inc -- ускорение) или dv_X = inc (inc_is_velocity, для
// XSPH: inc = av); v_mid = base + (h*dvx + av)/2, av -- если upd.av задан.
inline double work_mid_solid(const rhs::CpuExecutor& exec, const ParticleSystem& sys,
                             const ode::SolidVelocityUpdate& upd, const std::vector<std::vector<double>>& dvx,
                             std::size_t n, const std::vector<std::vector<double>>& inc, bool inc_is_velocity)
{
    const double k = inc_is_velocity ? 1.0 : upd.h;
    return exec.parallel_reduce(
        n, 0.0,
        [&](std::size_t i) -> double {
            double w = 0.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                double dv_step = upd.h * dvx[d][i];
                if (upd.av != nullptr) dv_step += (*upd.av)[d][i];
                const double v_mid = (*upd.base[d])[i] + 0.5 * dv_step;
                w += k * inc[d][i] * v_mid;
            }
            return energy_mass(sys, i) * w;
        },
        [](double a, double b) { return a + b; });
}

// Sum_i energy_mass_i * x_i для i < n (в 2D осесимметрии -- масса кольца).
inline double mass_weighted_sum(const rhs::CpuExecutor& exec, const ParticleSystem& sys,
                                const std::vector<double>& x, std::size_t n)
{
    return exec.parallel_reduce(
        n, 0.0, [&](std::size_t i) -> double { return energy_mass(sys, i) * x[i]; },
        [](double a, double b) { return a + b; });
}

// CSV-журнал: одна строка на записываемый шаг.
class EnergyLogger
{
public:
    // Открывает CSV (перезапись), пишет заголовок; бросает std::runtime_error,
    // если файл не открылся.
    // append = true (рестарт из чекпоинта): файл дописывается, без заголовка;
    // базу E_mech0 вызывающий восстанавливает через restore_e0().
    EnergyLogger(const std::string& path, std::size_t every, bool append = false) : every_(every)
    {
        out_.open(path, append ? (std::ios::out | std::ios::app) : (std::ios::out | std::ios::trunc));
        if (!out_.is_open()) throw std::runtime_error("EnergyLogger: не удалось открыть файл " + path);
        out_ << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (append) return;
        out_ << "step,time,dt,KE_solid,KE_fiber,U_int,U_shear_el,U_bond,U_hertz_ff,U_hertz_fs,"
                "E_mech,D_friction,D_xsph,D_contact_ss,D_break,S_u_clamp,Q_av,W_plastic,residual,"
                "W_p_fiber,W_v_fiber,D_break_contact,clip_loss,D_contact_fiber\n";
        out_.flush();
    }

    std::size_t every() const { return every_; }
    bool due(std::size_t step) const { return every_ > 0 && step % every_ == 0; }

    EnergyStockWorkspace& workspace() { return ws_; }

    // База энергии для чекпоинта.
    double e_mech0() const { return e_mech0_; }
    bool have_e0() const { return have_e0_; }
    void restore_e0(double e_mech0, bool have_e0)
    {
        e_mech0_ = e_mech0;
        have_e0_ = have_e0;
    }

    // Ускорение трения шага: перед стадией трения сюда копируется dv/dvx, после
    // неё (и повторной проекции зажимов) -- заменяется разностью «после - до».
    // Память только растёт.
    std::vector<double> friction_fiber[dim];
    std::vector<std::vector<double>> friction_solid;
    bool has_friction = false; // friction_* заполнены на этом шаге

    // Ускорение контактов нитей с законом 1-3 (Fiber-Fiber и Fiber-Solid, сторона нити;
    // Fiber-Solid, сторона Solid) -- для канала D_contact_fiber. Память только растёт.
    std::vector<double> contact_fiber_acc[dim];
    std::vector<std::vector<double>> contact_solid_acc;
    bool has_contact_fiber = false; // contact_*_acc заполнены на этом шаге

    // Первая записанная строка задаёт E_mech0;
    // residual = (E_mech - E_mech0) + friction + xsph + contact_ss + brk - u_clamp
    //            + w_p_fiber + w_v_fiber + brk_contact - clip + contact_fiber.
    void write_row(std::size_t step, double time, double dt, const EnergyStocks& s, const EnergyChannels& c)
    {
        const double e_mech = s.e_mech();
        if (!have_e0_)
        {
            e_mech0_ = e_mech;
            have_e0_ = true;
        }
        const double residual = (e_mech - e_mech0_) + c.friction + c.xsph + c.contact_ss + c.brk - c.u_clamp +
                                c.w_p_fiber + c.w_v_fiber + c.brk_contact - c.clip + c.contact_fiber;
        out_ << step << ',' << time << ',' << dt << ',' << s.ke_solid << ',' << s.ke_fiber << ',' << s.u_int << ','
             << s.u_shear_el << ',' << s.u_bond << ',' << s.u_hertz_ff << ',' << s.u_hertz_fs << ',' << e_mech << ','
             << c.friction << ',' << c.xsph << ',' << c.contact_ss << ',' << c.brk << ',' << c.u_clamp << ','
             << c.q_av << ',' << c.w_plastic << ',' << residual << ',' << c.w_p_fiber << ',' << c.w_v_fiber << ',' << c.brk_contact << ',' << c.clip << ','
             << c.contact_fiber << '\n';
        out_.flush();
    }

private:
    std::ofstream out_;
    std::size_t every_ = 0;
    double e_mech0_ = 0.0;
    EnergyStockWorkspace ws_;
    bool have_e0_ = false;
};

} // namespace energy

#endif // SPH_ENERGY_LOG
