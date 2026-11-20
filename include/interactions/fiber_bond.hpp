#pragma once
// interactions/fiber_bond.hpp
// Fiber-Bond: сила упругой связи между топологическими соседями одной нити
// (notes/fiber-model/fiber-problem-statement.md, раздел «Fiber-Bond», «Постановка проекта: только
// растяжение»; notes/fiber-model/fiber-plan.md, этап Ф3, §3.2/§3.3/§3.7).
//
// Primitive — уже существующая fiber::FiberBond (i, j, l0, l_break),
// fiber/fiber_system.hpp, этап Ф2: список связей строится из topology-графа
// нити (FiberTopology::bonds), а не пространственным поиском (fiber-plan.md §3.2,
// таблица "Типы взаимодействий и primitives").
//
// Kernel только ЧИТАЕТ FiberState/FiberMaterial (readonly) и производит
// typed contributions (FiberDelta) через rhs::ContributionWriter —
// суммирование происходит отдельно, в rhs::reduce_contributions. Никаких
// atomics/mutex здесь нет (notes/fiber-model/fiber-problem-statement.md, «Основные требования к interaction
// kernel», пп. 1-5).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "fiber/fiber_plasticity.hpp"
#include "fiber/fiber_system.hpp"
#include "param.hpp"
#include "rhs/contributions.hpp"
#include "rhs/cpu_executor.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Метка interaction_type для ContributionMeta/детерминированного reduce
// (rhs/contributions.hpp, «Детерминизм») — произвольная константа, уникальная
// среди типов взаимодействий, вносящих вклад в Fiber.
constexpr std::uint32_t kFiberBondInteractionType = 1;

// Площадь поперечного сечения нити (notes/fiber-model/fiber-problem-statement.md, Fiber-Bond): A = pi*d^2/4.
// ВНИМАНИЕ: модуль Юнга нити — FiberMaterial.E, не FiberMaterial.nu
// (коэффициент Пуассона, в Fiber-Bond не используется — постановка проекта
// работает только с осевой жёсткостью EA, изгиб не моделируется).
inline double fiber_cross_section_area(const fiber::FiberMaterial& mat)
{
    return pi * mat.d * mat.d / 4.0;
}

// Эффективная мгновенная жёсткость E_inst = E + E1 (ветвь Максвелла не успела
// расслабиться; дизайн, раздел 5). Ветвь учитывается, только если задана
// полностью (E1 > 0 и tau > 0), как в fiber::bond_constitutive_update.
inline double fiber_instant_modulus(const fiber::FiberMaterial& mat)
{
    return mat.E + ((mat.visco.E1 > 0.0 && mat.visco.tau > 0.0) ? mat.visco.E1 : 0.0);
}

// Вязкость Кельвина ЭТОЙ связи. Если задана eta [Па·с] -- она и есть. Если задан
// xi (доля критического демпфирования, имеет приоритет над eta), eta считается
// через fiber::kelvin_eta_from_xi по E_inst и приведённой массе бусин связи.
// Считается на лету (одна sqrt и деление) и в RHS-ядре, и в проходе
// Postprocessing, а не кэшируется: масса бусин и l0 постоянны, значит результат
// один и тот же в обоих местах (детерминированно), а кэш потребовал бы ещё
// одного поля связи, которое надо переносить при разрыве. Возвращает копию
// ViscoParams с заполненным eta -- fiber::BondLaw хранит ссылки.
inline fiber::ViscoParams fiber_bond_visco(const fiber::FiberMaterial& mat, double A, double l0, double mi,
                                           double mj)
{
    fiber::ViscoParams vp = mat.visco;
    if (vp.xi > 0.0)
    {
        const double m_red = (mi * mj) / (mi + mj);
        vp.eta = fiber::kelvin_eta_from_xi(vp.xi, fiber_instant_modulus(mat), A, l0, m_red);
    }
    return vp;
}

// Участок разупрочнения Δε_f ЭТОЙ связи (notes/fiber-model/fiber-damage-plan.md): DEPS_F напрямую или
// из энергии разрушения GF: Δε_f = 2·GF/(σ_u·l0) -- энергия разрыва связи не зависит от l0.
// 0, если повреждённость выключена.
inline double fiber_bond_damage_deps_f(const fiber::FiberMaterial& mat, double l0)
{
    if (!mat.damage.enabled() || !(mat.failure.sigma_u > 0.0)) return 0.0;
    if (mat.damage.deps_f > 0.0) return mat.damage.deps_f;
    return l0 > 0.0 ? 2.0 * mat.damage.gf / (mat.failure.sigma_u * l0) : 0.0;
}

// Критический безразмерный шаг h = omega*dt схемы проекта (предиктор -> RHS ->
// корректор) для линейного осциллятора с демпфером Кельвина, доля критического
// демпфирования xi: h_cr = 2 (sqrt(1+4 xi^2) - 2 xi). xi = 0 даёт 2, при
// больших xi -> 1/(2 xi). Вторая форма без вычитания близких чисел.
inline double fiber_kelvin_h_crit(double xi)
{
    return 2.0 / (std::sqrt(1.0 + 4.0 * xi * xi) + 2.0 * xi);
}

// ----------------------------------------------------------------------
// compute_fiber_bond_contributions — заполняет thread_buffers (по одному
// на CPU worker, как в rhs/contributions.hpp) contributions Fiber-Bond для
// ВСЕХ связей topology.bonds. Буферы переиспользуются между вызовами
// (изменяют размер под exec.num_workers() при необходимости, иначе только
// clear() — «память только растёт», notes/fiber-model/fiber-problem-statement.md «Временная память»).
//
// Для каждой связи (i, j, l0, l_break):
//   x_ij = x_j - x_i,  l_ij = |x_ij|
//   F_i<-j = E*A * max(l_ij/l0 - 1, 0) * x_ij/l_ij      (только растяжение)
//   F_j<-i = -F_i<-j
//   dv_i = F_i<-j / m_i,  dv_j = F_j<-i / m_j = -F_i<-j / m_j
// Для материала с неупругостью (FiberMaterial::inelastic(), дизайн
// notes/fiber-model/fiber-plasticity-design.md, раздел 3) вместо E*strain берётся
// fiber::bond_stress: sigma = max(E*(eps - eps_p) + q + eta_v*eps_dot, 0),
// сила A*sigma вдоль связи, eps_dot = ((v_j - v_i).n)/l0 по текущим скоростям
// FiberState::vx (на момент RHS), eps_p и q читаются из FiberBond::st. Ядро
// состояние НЕ пишет. Сила остаётся центральной и противоположной -- импульс и
// момент сохраняются. Для упругого материала (по умолчанию) выполняется ровно
// прежняя ветка с прежним порядком операций (побитовая совместимость, Т0).
// l_break здесь не используется — критерий разрыва проверяется в
// Postprocessing (этап Ф7, notes/fiber-model/fiber-problem-statement.md «Разрыв связи»), НЕ в этом kernel.
// ----------------------------------------------------------------------
inline void compute_fiber_bond_contributions(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers
#ifdef SPH_ENERGY_LOG
    ,
    // Лог энергии: напряжение модели σ_m [Па] по связям там, где отсечка σ >= 0 сработала
    // при ε > ε_p (σ_m < 0), иначе 0; индекс = индекс связи. Только запись (по k, без гонок).
    std::vector<double>* sigma_clip = nullptr
#endif
)
{
    const std::size_t nworkers = exec.det_segments(); // детерминированно: не зависит от числа потоков
    if (thread_buffers.size() < nworkers) thread_buffers.resize(nworkers);
    for (auto& buf : thread_buffers) buf.clear();

    const std::size_t nbonds = topology.bonds.size();
#ifdef SPH_ENERGY_LOG
    if (sigma_clip) sigma_clip->assign(nbonds, 0.0);
#endif
    if (nbonds == 0) return;

    exec.parallel_for_segments(nbonds, [&](std::size_t seg, std::size_t k) {
        const fiber::FiberBond& b = topology.bonds[k];
        const std::size_t i = b.i;
        const std::size_t j = b.j;
        const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(fiber::bond_material_id(b, state))];
        const double A = fiber_cross_section_area(mat);

        double xij[dim];
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            xij[d] = state.x[d][j] - state.x[d][i];
            l2 += xij[d] * xij[d];
        }
        const double l = std::sqrt(l2);

        fiber::FiberDelta delta_i{};
        fiber::FiberDelta delta_j{};

        const double mi = state.mass[i];
        const double mj = state.mass[j];
        bool active = false;
        double coeff = 0.0; // F_i<-j = coeff * xij
        if (mat.inelastic())
        {
            if (l > 0.0)
            {
                const double strain = l / b.l0 - 1.0;
                double eps_dot = 0.0;
                const fiber::ViscoParams vp = fiber_bond_visco(mat, A, b.l0, mi, mj);
                if (vp.eta != 0.0)
                {
                    double vrel = 0.0; // (v_j - v_i) . x_ij
                    for (std::size_t d = 0; d < dim; ++d) vrel += (state.vx[d][j] - state.vx[d][i]) * xij[d];
                    eps_dot = (vrel / l) / b.l0;
                }
                const fiber::BondLaw law(mat.E, mat.plastic, mat.rate, vp, mat.failure, mat.damage, fiber_bond_damage_deps_f(mat, b.l0));
                const double sigma = fiber::bond_stress(law, b.st, strain, eps_dot);
#ifdef SPH_ENERGY_LOG
                if (sigma_clip && strain > b.st.eps_p)
                {
                    const double sigma_m = law.E * (strain - b.st.eps_p) + b.st.q + (vp.eta != 0.0 ? vp.eta * eps_dot : 0.0);
                    if (sigma_m < 0.0) (*sigma_clip)[k] = sigma_m;
                }
#endif
                coeff = A * sigma / l;
                active = true;
            }
        }
        else if (l > b.l0 && l > 0.0) // только растяжение (notes/fiber-model/fiber-problem-statement.md, "Постановка проекта")
        {
            const double strain = l / b.l0 - 1.0;
            coeff = mat.E * A * strain / l;
            active = true;
        }
        if (active)
        {
            for (std::size_t d = 0; d < dim; ++d)
            {
                const double f = coeff * xij[d];
                delta_i.dv[d] = f / mi;  // +F_i<-j / m_i
                delta_j.dv[d] = -f / mj; // -F_i<-j / m_j == F_j<-i / m_j
            }
        }

        rhs::ContributionWriter<fiber::FiberDelta> writer(thread_buffers[seg],
                                                            kFiberBondInteractionType);
        writer.emit(i, delta_i, /*primitive_id=*/static_cast<std::uint64_t>(k), /*slot=*/0);
        writer.emit(j, delta_j, /*primitive_id=*/static_cast<std::uint64_t>(k), /*slot=*/1);
    });
}

// compute_fiber_bond — удобная обёртка: contributions + reduce за один
// вызов. Записывает в output вектор FiberDelta длины state.size() (индекс —
// частица Fiber). thread_buffers, output и scratch передаются вызывающей
// стороной (обычно — переиспользуемые переменные цикла по времени) для
// соблюдения "память только растёт".
inline void compute_fiber_bond(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                               const fiber::FiberTopology& topology,
                               const std::vector<fiber::FiberMaterial>& materials,
                               std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers,
                               std::vector<fiber::FiberDelta>& output, rhs::ReduceScratch& scratch
#ifdef SPH_ENERGY_LOG
                               ,
                               std::vector<double>* sigma_clip = nullptr
#endif
)
{
#ifdef SPH_ENERGY_LOG
    compute_fiber_bond_contributions(exec, state, topology, materials, thread_buffers, sigma_clip);
#else
    compute_fiber_bond_contributions(exec, state, topology, materials, thread_buffers);
#endif
    rhs::reduce_contributions_into(exec, thread_buffers, state.size(), output, scratch);
}

// Обёртка со старой сигнатурой: результат по значению, локальная рабочая память.
inline std::vector<fiber::FiberDelta> compute_fiber_bond(
    const rhs::CpuExecutor& exec, const fiber::FiberState& state, const fiber::FiberTopology& topology,
    const std::vector<fiber::FiberMaterial>& materials,
    std::vector<rhs::ContributionBuffer<fiber::FiberDelta>>& thread_buffers)
{
    std::vector<fiber::FiberDelta> output;
    rhs::ReduceScratch scratch;
    compute_fiber_bond(exec, state, topology, materials, thread_buffers, output, scratch);
    return output;
}

// Член минимума compute_fiber_bond_dt_limit для одной связи (формула -- в шапке
// compute_fiber_bond_dt_limit ниже: упругий 0.8/omega_inst и, при демпфере
// Кельвина, предел Белычко по старшей моде цепочки). Общий с dt_ref нити в
// законах контакта (interactions/fiber_contact_law.hpp, compute_fiber_bead_dt_ref).
inline double fiber_bond_pair_dt(const fiber::FiberMaterial& mat, double l0, double mi, double mj)
{
    const double A = fiber_cross_section_area(mat);
    const double k_eff = fiber_instant_modulus(mat) * A / l0; // E1 = 0 => ровно прежнее E*A/l0
    const double m_red = (mi * mj) / (mi + mj);
    const double omega = std::sqrt(k_eff / m_red);
    double dt_b = 0.8 / omega;
    if (mat.visco.xi > 0.0 || mat.visco.eta > 0.0)
    {
        const double xi = mat.visco.xi > 0.0 ? mat.visco.xi
                                             : (mat.visco.eta * A / l0) / (2.0 * std::sqrt(k_eff * m_red));
        // Старшая мода цепочки: omega*sqrt2, xi*sqrt2 (см. шапку compute_fiber_bond_dt_limit).
        constexpr double kSqrt2 = 1.4142135623730951;
        dt_b = std::min(dt_b, 0.8 * fiber_kelvin_h_crit(kSqrt2 * xi) / (kSqrt2 * omega));
    }
    return dt_b;
}

// ----------------------------------------------------------------------
// compute_fiber_bond_dt_limit — dt_limit связей (notes/fiber-model/fiber-problem-statement.md, "Шаг по
// времени"; notes/fiber-model/fiber-plan.md §3.7):
//     dt_lim = 0.8 / sqrt(k_eff / m_red),  k_eff = E*A/l0,
//     m_red = m_i*m_j / (m_i+m_j).
// Для неупругого материала (дизайн, раздел 5): жёсткость по E_inst = E + E1,
// а при демпфере Кельвина множитель Белычко (sqrt(1+xi^2) - xi),
//     dt_lim = min( 0.8/omega_inst,  0.8 * h_cr(sqrt2*xi) / (sqrt2*omega_inst) ),
//     h_cr(x) = 2 (sqrt(1+4x^2) - 2x) = 2/(sqrt(1+4x^2) + 2x),
//     xi = c/(2 sqrt(k_inst m_red)), c = eta_v*A/l0 (или xi из материала).
// h_cr(xi) -- точный критический шаг h = omega*dt схемы проекта на паре
// (предиктор скорости -> RHS с демпфером по скорости предиктора -> корректор;
// спектральный радиус матрицы шага 3x3 выходит на 1 через -1: это формула
// Белычко с удвоенным ξ из-за запаздывающей скорости предиктора). Старшая мода
// ЦЕПОЧКИ (бусина в двух связях, оба конца связи по m_red = m/2) имеет
// omega_max = sqrt2*omega_пары и, так как c ~ k, xi_max = sqrt2*xi_пары -- шаг
// считается по ней, а не по паре (иначе при xi >~ 1.2 цепочка неустойчива даже
// при 0.9*dt_lim). Множитель 0.8 -- тот же запас, что в упругом лимите, а min
// с упругим 0.8/omega_inst гарантирует, что демпфер лимит никогда не ослабляет.
// Подтверждено тестом Т9 (пара и цепочка, tests/test_fiber_plasticity_integration.cpp).
// Без демпфера Кельвина (в т.ч. при одной пластике) значение прежнее побитово.
// Считается для ВСЕХ связей независимо от знака деформации ("чтобы шаг не
// прыгал при переходе через l0") — минимум по всем связям через
// exec.parallel_reduce. Возвращает +inf (numeric_limits::max), если связей
// нет — так вызывающая сторона может безусловно брать min(..., это
// значение) без отдельной проверки на пустоту.
// ----------------------------------------------------------------------
inline double compute_fiber_bond_dt_limit(const rhs::CpuExecutor& exec, const fiber::FiberState& state,
                                           const fiber::FiberTopology& topology,
                                           const std::vector<fiber::FiberMaterial>& materials)
{
    const std::size_t nbonds = topology.bonds.size();
    const double kInf = std::numeric_limits<double>::max();
    if (nbonds == 0) return kInf;

    return exec.parallel_reduce(
        nbonds, kInf,
        [&](std::size_t k) -> double {
            const fiber::FiberBond& b = topology.bonds[k];
            const fiber::FiberMaterial& mat = materials[static_cast<std::size_t>(fiber::bond_material_id(b, state))];
            return fiber_bond_pair_dt(mat, b.l0, state.mass[b.i], state.mass[b.j]);
        },
        [](double a, double b2) { return std::min(a, b2); });
}

} // namespace interactions
