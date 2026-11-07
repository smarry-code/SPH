#pragma once
// guo2019_monitor.hpp
// Монитор шага и критерий исхода выстрела для опытов Guo 2019 (ТЗ Guo 2019 §4.1,
// §6.1-§6.3). GuoShotMonitor -- наблюдатель RunGuard::set_step_observer: после
// каждого шага читает fiber::Model (по const-ссылке, состояние НЕ меняет),
// пишет monitor.csv (§6.1), определяет исход (§6.2) и по исходу просит у сторожа
// штатный останов (request_stop). Сводка выстрела -- write_shot_txt (§6.3).
//
// От ctest не зависит (нет e2e::check): им пользуются и smoke-тест
// tests/test_e2e_guo2019.cpp, и программа выстрелов tools/guo2019/guo2019_shot.cpp.
//
// Что считается на каждом шаге (O(N)): энергии, z_c/v_z шара, зажим (S6), суммы
// масс (S7), число связей, прогибы. Дорогое (контакты нить--нить, нить--Solid,
// туннелирование, m_ej) -- раз в n_diag шагов, а также на шаге 1 и в finalize().
//
// Энергия. E_tot = KE_s + KE_p + KE_f + U_b + U_int, где U_int = sum m u (Solid)
// минус стартовое значение. Члена работы зажима нет: у неподвижного зажима ядра
// она равна нулю тождественно (ТЗ Guo 2019 §4.1). Энергия Герца в контактах и
// энергия порванных связей в E_tot не входят (известная неполнота, §6.1).

#include "guo2019_scene.hpp"

#include "clamp_e2e_common.hpp"
#include "fiber/model.hpp"
#include "interactions/fiber_contact.hpp"
#include "interactions/solid_contact.hpp"
#include "rhs/cpu_executor.hpp"
#include "run_guard.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h> // fsync

namespace guo2019
{

// Штатный останов по запросу монитора (ТЗ Guo 2019 §4.1). Наблюдатель вызывается
// в RunGuard::should_stop() первым делом, а счётчик шагов steps_ увеличивается сразу
// после него. Поэтому set_step_limit(steps() + 1) делает так, что should_stop()
// вернёт true с причиной StepLimit на ЭТОМ ЖЕ шаге. Аварией это не считается
// (failed() == false). Если на том же шаге сработал WallClock или Stalled, причиной
// останется авария -- так и задумано (StepLimit проверяется последним).
inline void request_stop(RunGuard& g)
{
    g.set_step_limit(g.steps() + 1);
}

enum class Outcome
{
    Running,       // исход ещё не определён
    Perforated,    // пробит: свободный полёт за мишенью с v_z < 0
    NotPerforated, // не пробит: v_z >= 0 непрерывно kOutcomeWindow
    Undetermined   // достигнут t_max (или расчёт закончен) без исхода
};

inline const char* outcome_name(Outcome o)
{
    switch (o)
    {
    case Outcome::Running: return "undetermined";
    case Outcome::Perforated: return "perforated";
    case Outcome::NotPerforated: return "not_perforated";
    case Outcome::Undetermined: return "undetermined";
    }
    return "?";
}

struct MonitorOptions
{
    std::size_t n_mon = 1;       // строка monitor.csv раз в n_mon шагов
    std::size_t n_diag = 1;      // дорогая диагностика раз в n_diag шагов
    double t_max = 1.0e30;       // потолок времени до исхода, с
    std::string csv_path;        // пусто -- csv не пишется
    bool stop_on_outcome = true; // останавливать расчёт при определённом исходе или t >= t_max
    RunGuard* guard = nullptr;   // для request_stop
    double wp_radius = 0.0;      // круг для w_p; <= 0 -- sc.spec.wp_radius (30 мм, ТЗ Guo 2019 §6.1; scaled 7.5 мм)
    // Продолжение из чекпоинта: monitor.csv не перезаписывается в конструкторе; restore_state() оставляет в нём
    // строки до шага чекпоинта и дописывает дальше (без restore_state -- файл заводится заново на первом шаге).
    bool restart = false;
};

// Все накопленные величины (читаются после finalize()).
struct MonitorStats
{
    std::size_t steps = 0;
    double t_end = 0.0;
    double dt_first = 0.0, dt_min = std::numeric_limits<double>::infinity();

    // Энергии
    double E_tot_0 = 0.0, E_tot_max = 0.0, E_tot_end = 0.0;
    double KE_s_0 = 0.0, KE_s_end = 0.0, KE_p_end = 0.0, KE_f_end = 0.0, U_b_end = 0.0, U_int_end = 0.0;

    // Шар
    double z_c_end = 0.0, v_z_end = 0.0;

    // Реакции зажимов: наибольшая по модулю z-компонента (со знаком) и момент
    double F_s_z_max = 0.0, F_s_z_max_t = 0.0;
    double F_f_z_max = 0.0, F_f_z_max_t = 0.0;
    std::size_t n_clamped_s = 0, n_clamped_f = 0;

    // Нити
    std::size_t n_bonds_0 = 0, n_bonds_end = 0;
    std::size_t n_broken_v = 0, n_broken_h = 0;
    double t_first_break = -1.0;       // < 0 -- разрывов не было
    double w_f_at_first_break = 0.0;
    double w_f_max = 0.0, w_p_max = 0.0;
    double w_f_end = 0.0;
    // Неупругая нить (include/fiber/fiber_plasticity.hpp): максимум eps_p по связям и шагам, суммарная пластическая
    // работа W_p [Дж] по всем связям, включая порванные (значение связи на последнем шаге перед удалением).
    double eps_p_max = 0.0, W_p_total = 0.0;
    std::size_t n_bonds_plastic = 0; // связей с eps_p > 0 на последнем шаге

    // Пакет слоёв (n_plies >= 2; ТЗ 2T/4T §4.3). При n_plies == 1 не заполняется и не пишется.
    int n_plies = 1;
    std::vector<std::size_t> n_broken_v_p, n_broken_h_p;      // порванные нити по слоям (слой 1 -- лицевой)
    std::vector<double> t_first_break_p;                      // < 0 -- у слоя разрывов не было
    std::vector<double> w_f_at_first_break_p;                 // прогиб слоя в момент его первого разрыва
    std::vector<double> pullin_at_first_break_p;              // втягивание свободных кромок слоя при его первом разрыве
    std::vector<double> pullin_end_p;                         // то же на t_end (оценка распускания, не измерение)
    std::vector<double> w_f_end_p;                            // прогиб слоя на последнем шаге
    double contacts_per_crossing_inter = 0.0;                 // максимум по шагам, на (N - 1) nx ny
    double contacts_per_crossing_inter_step1 = 0.0;
    double max_delta_ff_inter_d = 0.0;                        // максимум по шагам, контакты между слоями
    std::size_t ff_intra_initial = 0, ff_inter_initial = 0;   // контакты нить--нить на исходной укладке
    std::size_t ff_intra_step1 = 0, ff_inter_step1 = 0;       // на шаге 1
    std::size_t ff_intra_last = 0, ff_inter_last = 0;         // на последнем шаге дорогой диагностики
    std::size_t ply_order_violations = 0;                     // максимум по шагам диагностики
    std::size_t ply_order_violations_step1 = 0;
    double t_ply_order_first = -1.0;                          // < 0 -- нарушений не было
    double ply_gap_min = 0.0, ply_gap_max = 0.0;              // зазор слоёв по вершинным бусинам на старте (ply_vertex_gaps)

    // Дорогая диагностика
    double m_ej = 0.0, m_ej_max = 0.0;
    std::size_t n_tunnel = 0;                 // максимум по шагам
    double max_delta_ff_d = 0.0;              // максимум по шагам
    double max_delta_fs_d = 0.0;
    double contacts_per_crossing = 0.0;       // максимум по шагам
    double contacts_per_crossing_step1 = 0.0; // на шаге 1
    std::size_t ff_contacts_initial = 0;      // принятые контакты Fiber-Fiber на исходной укладке (до шага 1)
    double ff_delta_d_initial = 0.0;
    std::size_t ff_contacts_step1 = 0;
    double ff_delta_d_step1 = 0.0;
    std::size_t ff_contacts_max = 0;

    // S2: первый шаг
    double plate_vmax_step1 = 0.0;
    double bead_vmax_step1 = 0.0;

    // S6 / S7
    std::size_t clamp_violations = 0;         // число (шаг, элемент) нарушений неподвижности
    std::string clamp_first_violation;
    std::size_t mass_violations = 0;          // число шагов, где сумма масс отличалась от стартовой

    // Исход
    Outcome outcome = Outcome::Running;
    double V_R = 0.0;
    double t_outcome = 0.0;
    bool stopped_by_monitor = false;          // монитор вызвал request_stop

    // Собственная стоимость монитора
    double monitor_wall_s = 0.0;

    // Проекция скорости нормального контакта (fiber_*_contact_predictor = 1): максимумы за ВЕСЬ
    // прогон по шагам (из fiber::Model::proj_stats). При выключенной проекции не заполняются.
    double proj_resid_max = 0.0, proj_resid_rel_max = 0.0;
    double W_proj_pos_max = 0.0;      // max(0, max W_stage): положительная работа стадии (теоретически 0)
    std::size_t proj_slots_max = 0;   // max (n_slots_ff + n_slots_fs) за шаг
};

namespace mon_detail
{

using clock = std::chrono::steady_clock;

inline double since(clock::time_point t0)
{
    return std::chrono::duration<double>(clock::now() - t0).count();
}

template <class T>
struct is_std_vector : std::false_type
{
};
template <class T, class A>
struct is_std_vector<std::vector<T, A>> : std::true_type
{
};

// Побайтовая (де)сериализация состояния монитора для блока USER чекпоинта: тривиально копируемые
// значения -- как есть, std::string и std::vector -- длина u64 и элементы.
struct StateOut
{
    std::string b;
    template <class T>
    void operator()(T& v)
    {
        if constexpr (std::is_same_v<T, std::string>)
        {
            std::uint64_t n = v.size();
            (*this)(n);
            b.append(v);
        }
        else if constexpr (is_std_vector<T>::value)
        {
            std::uint64_t n = v.size();
            (*this)(n);
            for (auto& e : v) (*this)(e);
        }
        else
        {
            static_assert(std::is_trivially_copyable_v<T>, "состояние монитора: неподдерживаемый тип");
            b.append(reinterpret_cast<const char*>(&v), sizeof(T));
        }
    }
};

struct StateIn
{
    const std::string& b;
    std::size_t pos = 0;
    void take(void* dst, std::size_t n)
    {
        if (n > b.size() - pos) throw std::runtime_error("состояние монитора в чекпоинте обрезано");
        std::memcpy(dst, b.data() + pos, n);
        pos += n;
    }
    template <class T>
    void operator()(T& v)
    {
        if constexpr (std::is_same_v<T, std::string>)
        {
            std::uint64_t n = 0;
            (*this)(n);
            if (n > b.size() - pos) throw std::runtime_error("состояние монитора в чекпоинте обрезано");
            v.assign(b.data() + pos, static_cast<std::size_t>(n));
            pos += static_cast<std::size_t>(n);
        }
        else if constexpr (is_std_vector<T>::value)
        {
            std::uint64_t n = 0;
            (*this)(n);
            if (n > b.size() - pos) throw std::runtime_error("состояние монитора в чекпоинте обрезано");
            v.resize(static_cast<std::size_t>(n));
            for (auto& e : v) (*this)(e);
        }
        else
        {
            static_assert(std::is_trivially_copyable_v<T>, "состояние монитора: неподдерживаемый тип");
            take(&v, sizeof(T));
        }
    }
};

} // namespace mon_detail

class GuoShotMonitor
{
public:
    GuoShotMonitor(const GuoScene& sc, const MonitorOptions& opt)
        : sc_(sc), m_(sc.model), opt_(opt), nt_(sc.ntotal), nb_(sc.model.fiber.size()),
          proj_on_(sc.cfg.fiber_fiber_contact_predictor == 1 || sc.cfg.fiber_solid_contact_predictor == 1)
    {
        if (!(opt_.wp_radius > 0.0)) opt_.wp_radius = sc.spec.wp_radius;
        const ParticleSystem& s = m_.solid;
        const fiber::FiberState& fb = m_.fiber;

        // Зажатые элементы: x0 (побитово).
        for (std::size_t i = 0; i < nt_ && i < m_.clamps.solid_mask.size(); ++i)
            if (m_.clamps.solid_mask[i] != 0) add_clamped(false, i, m_.clamps.solid_mask[i]);
        for (std::size_t k = 0; k < nb_ && k < m_.clamps.fiber_mask.size(); ++k)
            if (m_.clamps.fiber_mask[k] != 0) add_clamped(true, k, m_.clamps.fiber_mask[k]);

        // Массы
        mass_s0_ = sum_mass_solid();
        mass_f0_ = sum_mass_fiber();
        ball_mass_ = 0.0;
        for (std::size_t i = sc_.ball_begin; i < sc_.ball_end; ++i) ball_mass_ += s.mass[i];

        // Начальные z
        z0_bead_.resize(nb_);
        for (std::size_t k = 0; k < nb_; ++k) z0_bead_[k] = fb.x[2][k];
        z0_back_.resize(sc_.plate_back_layer.size());
        for (std::size_t q = 0; q < sc_.plate_back_layer.size(); ++q) z0_back_[q] = s.x[2][sc_.plate_back_layer[q]];

        // Нити: bead -> нить, число связей по нити в исходной топологии
        yarn_of_.assign(nb_, -1);
        for (std::size_t y = 0; y < sc_.yarns.size(); ++y)
            for (std::size_t k = sc_.yarns[y].bead_begin; k < sc_.yarns[y].bead_end && k < nb_; ++k)
                yarn_of_[k] = static_cast<int>(y);
        orig_bonds_.assign(sc_.yarns.size(), 0);
        for (std::size_t y = 0; y < sc_.yarns.size(); ++y)
            orig_bonds_[y] = sc_.yarns[y].bead_end - sc_.yarns[y].bead_begin - 1;
        cnt_.assign(sc_.yarns.size(), 0);
        st_.n_bonds_0 = m_.topology.bonds.size();
        last_bonds_ = st_.n_bonds_0;

        U_int0_ = clamp_e2e::internal_solid(s, nt_);
        const Energies e = energies();
        st_.KE_s_0 = e.ke_s;
        st_.E_tot_0 = e.ke_s + e.ke_p + e.ke_f + e.u_b; // U_int = 0 по определению
        st_.E_tot_max = st_.E_tot_0;
        st_.E_tot_end = st_.E_tot_0;

        // Пакет слоёв: бусина -> слой, вершинные бусины пересечений, концы нитей (ТЗ 2T/4T §4.3)
        init_plies();

        // Диаметр нити и число пересечений
        for (const auto& mat : m_.materials) d_fiber_ = std::max(d_fiber_, mat.d);
        crossings_ = static_cast<double>(sc_.has_fabric ? sc_.spec.nx * sc_.spec.ny : 0);

        // Fiber-Fiber на исходной укладке (до шага 1)
        if (nb_ > 0)
        {
            const ContactFF c = contact_ff();
            st_.ff_contacts_initial = c.count;
            st_.ff_delta_d_initial = c.delta_d;
            st_.ff_intra_initial = c.intra;
            st_.ff_inter_initial = c.inter;
        }

        if (!opt_.csv_path.empty())
        {
            if (opt_.restart)
                csv_pending_ = true; // открывается в restore_state() или на первом шаге
            else
                open_csv_fresh();
        }
    }

    // Новый monitor.csv с заголовком.
    void open_csv_fresh()
    {
        csv_pending_ = false;
        {
            csv_ = std::fopen(opt_.csv_path.c_str(), "w");
            if (csv_)
                std::fprintf(csv_,
                             "step,t,dt,z_c,v_z,KE_s,KE_p,KE_f,U_b,U_int,E_tot,Fs_x,Fs_y,Fs_z,Ff_x,Ff_y,Ff_z,"
                             "n_bonds,n_broken_v,n_broken_h,w_f,w_p");
            if (csv_)
            {
                if (nply_ >= 2)
                {
                    for (int k = 1; k <= nply_; ++k) std::fprintf(csv_, ",n_broken_v_p%d", k);
                    for (int k = 1; k <= nply_; ++k) std::fprintf(csv_, ",n_broken_h_p%d", k);
                    for (int k = 1; k <= nply_; ++k) std::fprintf(csv_, ",w_f_p%d", k);
                }
                // Только в режиме проекции: столбцы после всех прежних (при выключенной -- заголовок прежний).
                if (proj_on_)
                    std::fprintf(csv_, ",proj_n_ff,proj_n_fs,proj_max_ck,proj_max_cs,proj_resid,proj_resid_rel,"
                                       "proj_warm_frac,W_proj_pos_max");
                std::fprintf(csv_, "\n");
            }
        }
    }

    ~GuoShotMonitor()
    {
        if (csv_) std::fclose(csv_);
    }
    GuoShotMonitor(const GuoShotMonitor&) = delete;
    GuoShotMonitor& operator=(const GuoShotMonitor&) = delete;

    // Наблюдатель шага (RunGuard::set_step_observer).
    void on_step(double dt, double time)
    {
        const auto t0 = mon_detail::clock::now();
        const ParticleSystem& s = m_.solid;
        const fiber::FiberState& fb = m_.fiber;
        ++st_.steps;
        const std::size_t step = st_.steps;
        if (csv_pending_)
        {
            std::fprintf(stderr, "ПРЕДУПРЕЖДЕНИЕ: монитор: рестарт без состояния монитора -- monitor.csv заводится заново, "
                                 "вердикт и максимумы по шару только с шага рестарта.\n");
            open_csv_fresh();
        }
        st_.t_end = time;
        if (st_.dt_first <= 0.0 && dt > 0.0) st_.dt_first = dt;
        if (dt > 0.0) st_.dt_min = std::min(st_.dt_min, dt);

        Row r;
        r.step = step;
        r.t = time;
        r.dt = dt;

        // --- шар: центр масс и скорость
        {
            double zm = 0.0, vm = 0.0;
            for (std::size_t i = sc_.ball_begin; i < sc_.ball_end; ++i)
            {
                zm += s.mass[i] * s.x[2][i];
                vm += s.mass[i] * s.vx[2][i];
            }
            r.z_c = zm / ball_mass_;
            r.v_z = vm / ball_mass_;
        }

        // --- энергии
        const Energies e = energies();
        r.KE_s = e.ke_s;
        r.KE_p = e.ke_p;
        r.KE_f = e.ke_f;
        r.U_b = e.u_b;
        r.U_int = clamp_e2e::internal_solid(s, nt_) - U_int0_;
        r.E_tot = r.KE_s + r.KE_p + r.KE_f + r.U_b + r.U_int;
        st_.E_tot_max = std::max(st_.E_tot_max, r.E_tot);

        // --- реакции зажимов
        for (int d = 0; d < 3; ++d)
        {
            r.Fs[d] = m_.clamp_report.solid_force[d];
            r.Ff[d] = m_.clamp_report.fiber_force[d];
        }
        if (std::fabs(r.Fs[2]) > std::fabs(st_.F_s_z_max))
        {
            st_.F_s_z_max = r.Fs[2];
            st_.F_s_z_max_t = time;
        }
        if (std::fabs(r.Ff[2]) > std::fabs(st_.F_f_z_max))
        {
            st_.F_f_z_max = r.Ff[2];
            st_.F_f_z_max_t = time;
        }
        st_.n_clamped_s = m_.clamp_report.solid_clamped;
        st_.n_clamped_f = m_.clamp_report.fiber_clamped;

        // --- связи и разрывы
        r.n_bonds = m_.topology.bonds.size();
        if (r.n_bonds != last_bonds_)
        {
            last_bonds_ = r.n_bonds;
            count_broken(broken_v_, broken_h_);
        }
        r.n_broken_v = broken_v_;
        r.n_broken_h = broken_h_;
        if (nply_ >= 2)
        {
            r.pv = brk_pv_;
            r.ph = brk_ph_;
        }

        // --- пластическая работа и eps_p (только для неупругих материалов)
        accumulate_plastic();

        // --- прогибы
        {
            double wf = 0.0;
            for (std::size_t k = 0; k < nb_; ++k) wf = std::max(wf, z0_bead_[k] - fb.x[2][k]);
            r.w_f = nb_ > 0 ? wf : 0.0;
            double wp = 0.0;
            const double R2 = opt_.wp_radius * opt_.wp_radius;
            for (std::size_t q = 0; q < sc_.plate_back_layer.size(); ++q)
            {
                const std::size_t i = sc_.plate_back_layer[q];
                if (s.x[0][i] * s.x[0][i] + s.x[1][i] * s.x[1][i] <= R2) wp = std::max(wp, z0_back_[q] - s.x[2][i]);
            }
            r.w_p = wp;
        }
        st_.w_f_max = std::max(st_.w_f_max, r.w_f);
        st_.w_p_max = std::max(st_.w_p_max, r.w_p);
        if (st_.t_first_break < 0.0 && r.n_bonds < st_.n_bonds_0)
        {
            st_.t_first_break = time;
            st_.w_f_at_first_break = r.w_f;
        }
        if (nply_ >= 2) plies_step(r, time);

        // --- S2: скорости на первом шаге
        if (step == 1)
        {
            for (std::size_t i = sc_.plate_begin; i < sc_.plate_end; ++i)
                for (int d = 0; d < 3; ++d) st_.plate_vmax_step1 = std::max(st_.plate_vmax_step1, std::fabs(s.vx[d][i]));
            for (std::size_t k = 0; k < nb_; ++k)
                for (int d = 0; d < 3; ++d) st_.bead_vmax_step1 = std::max(st_.bead_vmax_step1, std::fabs(fb.vx[d][k]));
        }

        // --- S6: зажатые x == x0 побитово, v == 0.0
        for (const Clamped& c : clamped_)
        {
            for (int d = 0; d < 3; ++d)
            {
                if (!(c.mask & (1u << d))) continue;
                const double x = c.fiber ? fb.x[d][c.i] : s.x[d][c.i];
                const double v = c.fiber ? fb.vx[d][c.i] : s.vx[d][c.i];
                if (!(x == c.x0[d]) || !(v == 0.0))
                {
                    ++st_.clamp_violations;
                    if (st_.clamp_first_violation.empty())
                    {
                        char buf[256];
                        std::snprintf(buf, sizeof buf, "шаг %zu, %s %zu, компонента %d: x = %.17g (x0 = %.17g), v = %.17g",
                                      step, c.fiber ? "бусина" : "частица Solid", c.i, d, x, c.x0[d], v);
                        st_.clamp_first_violation = buf;
                    }
                }
            }
        }

        // --- S7: суммы масс побитово
        if (!(sum_mass_solid() == mass_s0_) || !(sum_mass_fiber() == mass_f0_)) ++st_.mass_violations;

        // --- дорогая диагностика
        if (step == 1 || (opt_.n_diag > 0 && step % opt_.n_diag == 0)) diagnostics(step == 1);

        // --- итоговые «последние» значения
        st_.E_tot_end = r.E_tot;
        st_.KE_s_end = r.KE_s;
        st_.KE_p_end = r.KE_p;
        st_.KE_f_end = r.KE_f;
        st_.U_b_end = r.U_b;
        st_.U_int_end = r.U_int;
        st_.z_c_end = r.z_c;
        st_.v_z_end = r.v_z;
        st_.w_f_end = r.w_f;
        st_.n_bonds_end = r.n_bonds;
        st_.n_broken_v = r.n_broken_v;
        st_.n_broken_h = r.n_broken_h;
        if (nply_ >= 2)
        {
            st_.n_broken_v_p.assign(r.pv.begin(), r.pv.end());
            st_.n_broken_h_p.assign(r.ph.begin(), r.ph.end());
            st_.w_f_end_p = r.wp;
            st_.pullin_end_p = pullin_cur_;
        }

        // --- исход
        update_outcome(time, r.v_z, r.z_c);

        // --- проекция нормального контакта: накопление за окно между строками csv
        if (proj_on_ && m_.proj_stats.enabled) accumulate_proj(m_.proj_stats);

        // --- csv
        last_ = r;
        if (opt_.n_mon <= 1 || step == 1 || step % opt_.n_mon == 0)
        {
            write_row(r);
            last_written_ = true;
        }
        else
            last_written_ = false;

        st_.monitor_wall_s += mon_detail::since(t0);
    }

    // Вызывать после прогона: пересчитывает дорогую диагностику на конечном
    // состоянии (m_ej на t_end), дописывает последнюю строку csv, закрывает
    // исход (Running -> Undetermined).
    void finalize()
    {
        const auto t0 = mon_detail::clock::now();
        if (st_.steps > 0) diagnostics(false);
        if (st_.steps > 0 && !last_written_) write_row(last_);
        last_written_ = true;
        if (csv_)
        {
            std::fclose(csv_);
            csv_ = nullptr;
        }
        if (st_.outcome == Outcome::Running)
        {
            st_.outcome = Outcome::Undetermined;
            st_.t_outcome = st_.t_end;
        }
        st_.monitor_wall_s += mon_detail::since(t0);
    }

    const MonitorStats& stats() const { return st_; }

    // ---- состояние для чекпоинта (блок USER, checkpoint::UserState) ----------------------------------
    // Зовётся при записи чекпоинта ПОСЛЕ on_step того же шага. Заодно сбрасывает monitor.csv на диск,
    // чтобы строки до шага чекпоинта пережили аварию.
    std::string save_state()
    {
        if (csv_)
        {
            std::fflush(csv_);
            fsync(fileno(csv_));
        }
        mon_detail::StateOut o;
        std::uint64_t magic = kStateMagic, nt = nt_, nb = nb_;
        std::int32_t nply = nply_;
        std::uint8_t proj = proj_on_ ? 1 : 0;
        o(magic);
        o(nt);
        o(nb);
        o(nply);
        o(proj);
        state_io(o);
        return std::move(o.b);
    }

    // После чтения чекпоинта (до первого шага). blob пуст -- состояния нет (старый чекпоинт): monitor.csv
    // заводится заново. Иначе состояние восстанавливается побитово, а в monitor.csv остаются строки до шага
    // чекпоинта включительно (после него -- дописываются). Несовпадение сцены или шага -- исключение.
    void restore_state(const std::string& blob, std::size_t guard_steps)
    {
        if (blob.empty())
        {
            if (csv_pending_)
            {
                std::fprintf(stderr, "ПРЕДУПРЕЖДЕНИЕ: в чекпоинте нет состояния монитора (формат v1/v2): monitor.csv "
                                     "заводится заново, вердикт и максимумы по шару только с шага рестарта.\n");
                open_csv_fresh();
            }
            return;
        }
        mon_detail::StateIn in{blob};
        std::uint64_t magic = 0, nt = 0, nb = 0;
        std::int32_t nply = 0;
        std::uint8_t proj = 0;
        in(magic);
        in(nt);
        in(nb);
        in(nply);
        in(proj);
        if (magic != kStateMagic) throw std::runtime_error("состояние монитора в чекпоинте: неизвестный формат");
        if (nt != nt_ || nb != nb_ || nply != nply_ || (proj != 0) != proj_on_)
            throw std::runtime_error("состояние монитора в чекпоинте от другой сцены (частицы, бусины, слои или проекция)");
        state_io(in);
        if (in.pos != blob.size()) throw std::runtime_error("состояние монитора в чекпоинте: лишние байты");
        if (st_.steps != guard_steps)
            throw std::runtime_error("состояние монитора в чекпоинте: шагов " + std::to_string(st_.steps) +
                                     ", у сторожа " + std::to_string(guard_steps));
        restored_ = true;
        if (csv_pending_) reopen_csv_after_restart();
    }
    bool state_restored() const { return restored_; }
    double ball_mass() const { return ball_mass_; }
    double ball_Reff() const { return sc_.R_eff; }
    double ball_z_c0() const { return sc_.z_c0; }
    // Энергия выстрела e = m_s (V0^2 - V_R^2)/2, ур. 1 статьи; NaN, если исход не определён.
    double shot_energy() const
    {
        if (st_.outcome != Outcome::Perforated && st_.outcome != Outcome::NotPerforated)
            return std::numeric_limits<double>::quiet_NaN();
        return 0.5 * ball_mass_ * (sc_.V0 * sc_.V0 - st_.V_R * st_.V_R);
    }

    // S15 (ТЗ 2T/4T §6.1): на ТЕКУЩЕМ состоянии считает контакты нить--нить двумя способами по одному
    // поиску: compute_fiber_contact_diagnostics (старый) и послойный обход примитивов (новый).
    struct FFSplitCheck
    {
        std::size_t old_count = 0, intra = 0, inter = 0;
        double old_delta_d = 0.0, delta_intra_d = 0.0, delta_inter_d = 0.0;
        bool match() const
        {
            return old_count == intra + inter && old_delta_d == std::max(delta_intra_d, delta_inter_d);
        }
    };
    FFSplitCheck ff_split_check() const
    {
        FFSplitCheck r;
        if (nb_ == 0 || nply_ < 2) return r;
        interactions::build_fiber_contact_search(diag_executor_, m_.fiber, m_.topology, m_.materials, diag_search_);
        const auto dg = interactions::compute_fiber_contact_diagnostics(diag_executor_, m_.fiber, m_.materials, diag_search_);
        r.old_count = dg.contact_count;
        r.old_delta_d = d_fiber_ > 0.0 ? dg.max_penetration / d_fiber_ : 0.0;
        const ContactFF c = split_contacts();
        r.intra = c.intra;
        r.inter = c.inter;
        r.delta_intra_d = c.delta_intra;
        r.delta_inter_d = c.delta_inter;
        return r;
    }

    // Разложение энергии (ТЗ Guo 2019 §6.2): пластина, ткань, шар, невязка.
    void print_energy_decomposition(std::FILE* f) const
    {
        const double E0 = st_.E_tot_0;
        std::fprintf(f, "Разложение энергии на t_end = %.6e с (E_tot_0 = %.6e Дж):\n", st_.t_end, E0);
        std::fprintf(f, "  шар (KE_s):                       %.6e Дж (%.2f %% от E_tot_0)\n", st_.KE_s_end,
                     100.0 * st_.KE_s_end / E0);
        std::fprintf(f, "  пластина (KE_p + U_int):          %.6e Дж (%.2f %%)\n", st_.KE_p_end + st_.U_int_end,
                     100.0 * (st_.KE_p_end + st_.U_int_end) / E0);
        std::fprintf(f, "  ткань (KE_f + U_b, без разрывов): %.6e Дж (%.2f %%)\n", st_.KE_f_end + st_.U_b_end,
                     100.0 * (st_.KE_f_end + st_.U_b_end) / E0);
        std::fprintf(f, "  невязка E_tot(t_end) - E_tot(0):  %.6e Дж (%.3f %%)\n", st_.E_tot_end - E0,
                     100.0 * (st_.E_tot_end - E0) / E0);
    }

    // shot.txt (ТЗ Guo 2019 §6.3): одна строка ключ=значение через пробел.
    void write_shot_txt(const std::string& path, double wall_s, const std::string& guard_reason) const
    {
        std::FILE* f = std::fopen(path.c_str(), "w");
        if (!f) throw std::runtime_error("guo2019: не удалось открыть " + path);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const bool det = st_.outcome == Outcome::Perforated || st_.outcome == Outcome::NotPerforated;
        const double VR = det ? st_.V_R : nan;
        const double E0 = st_.E_tot_0;
        std::fprintf(f, "config=%s level=%s V0=%.6g mu=%.6g outcome=%s V_R=%.9g e=%.9g t_end=%.9g steps=%zu wall_s=%.6g "
                        "dt_first=%.9g dt_min=%.9g E_tot_max/E_tot_0=%.9g E_tot_end/E_tot_0=%.9g "
                        "F_clamp_s_z_max=%.9g F_clamp_f_z_max=%.9g n_clamped_s=%zu n_clamped_f=%zu n_broken_v=%zu "
                        "n_broken_h=%zu w_f_max=%.9g t_first_break=%.9g w_f_at_first_break=%.9g w_p_max=%.9g "
                        "m_ej=%.9g n_tunnel=%zu max_delta_ff/d=%.9g max_delta_fs/d=%.9g contacts_per_crossing=%.9g "
                        "guard_reason=%s eps_p_max=%.9g W_p=%.9g n_bonds_plastic=%zu%s",
                     config_name(sc_.config), level_name(sc_.level), sc_.V0, sc_.mu, outcome_name(st_.outcome), VR,
                     shot_energy(), st_.t_end, st_.steps, wall_s, st_.dt_first, st_.dt_min, st_.E_tot_max / E0,
                     st_.E_tot_end / E0, st_.F_s_z_max, st_.F_f_z_max, st_.n_clamped_s, st_.n_clamped_f,
                     st_.n_broken_v, st_.n_broken_h, st_.w_f_max, st_.t_first_break < 0.0 ? nan : st_.t_first_break,
                     st_.t_first_break < 0.0 ? nan : st_.w_f_at_first_break, st_.w_p_max, st_.m_ej, st_.n_tunnel,
                     st_.max_delta_ff_d, st_.max_delta_fs_d, st_.contacts_per_crossing, guard_reason.c_str(), st_.eps_p_max,
                     st_.W_p_total, st_.n_bonds_plastic, (nply_ >= 2 || proj_on_) ? "" : "\n");
        if (nply_ >= 2)
        {
            // Новые ключи ТЗ 2T/4T §4.3 -- строго в конец строки (парсер tools/guo2019 читает по ключам).
            std::fprintf(f, " n_plies=%d", nply_);
            for (int k = 0; k < nply_; ++k)
            {
                const std::size_t q = static_cast<std::size_t>(k);
                const bool fb = st_.t_first_break_p[q] >= 0.0;
                std::fprintf(f,
                             " n_broken_v_p%d=%zu n_broken_h_p%d=%zu t_first_break_p%d=%.9g w_f_at_first_break_p%d=%.9g "
                             "pullin_at_first_break_p%d=%.9g pullin_end_p%d=%.9g",
                             k + 1, st_.n_broken_v_p[q], k + 1, st_.n_broken_h_p[q], k + 1, fb ? st_.t_first_break_p[q] : nan,
                             k + 1, fb ? st_.w_f_at_first_break_p[q] : nan, k + 1,
                             fb ? st_.pullin_at_first_break_p[q] : nan, k + 1, st_.pullin_end_p[q]);
            }
            std::fprintf(f,
                         " contacts_per_crossing_inter=%.9g max_delta_ff_inter/d=%.9g ply_order_violations=%zu "
                         "t_ply_order_first=%.9g ply_gap_min=%.12e ply_gap_max=%.12e%s",
                         st_.contacts_per_crossing_inter, st_.max_delta_ff_inter_d, st_.ply_order_violations,
                         st_.t_ply_order_first < 0.0 ? nan : st_.t_ply_order_first, st_.ply_gap_min, st_.ply_gap_max,
                         proj_on_ ? "" : "\n");
        }
        if (proj_on_)
        {
            // Ключи режима проекции -- строго в конец строки.
            std::fprintf(f, " proj_resid_max=%.9g proj_resid_rel_max=%.9g W_proj_pos_max=%.9g proj_slots_max=%zu\n",
                         st_.proj_resid_max, st_.proj_resid_rel_max, st_.W_proj_pos_max, st_.proj_slots_max);
        }
        std::fclose(f);
    }

private:
    // Накопление ε_p и W_p. Связь с w_p == 0 пропускается; W_p связи хранится по ключу (i,j) и после разрыва (удаления из
    // topology.bonds) остаётся в сумме со значением последнего шага.
    void accumulate_plastic()
    {
        bool any = false;
        for (const auto& mat : m_.materials)
            if (mat.inelastic()) { any = true; break; }
        if (!any) return;
        std::size_t npl = 0;
        for (const auto& b : m_.topology.bonds)
        {
            if (!(b.st.w_p > 0.0) && !(b.st.eps_p > 0.0)) continue;
            const fiber::FiberMaterial& mat = m_.materials[static_cast<std::size_t>(m_.fiber.material_id[b.i])];
            const double A = interactions::fiber_cross_section_area(mat);
            const double w = b.st.w_p * A * b.l0;
            const std::uint64_t key = (static_cast<std::uint64_t>(b.i) << 32) ^ static_cast<std::uint64_t>(b.j);
            double& slot = wp_bond_[key];
            st_.W_p_total += w - slot;
            slot = w;
            st_.eps_p_max = std::max(st_.eps_p_max, b.st.eps_p);
            if (b.st.eps_p > 0.0) ++npl;
        }
        st_.n_bonds_plastic = npl;
    }
    std::unordered_map<std::uint64_t, double> wp_bond_;

    struct Clamped
    {
        bool fiber;
        std::size_t i;
        std::uint8_t mask;
        double x0[3];
    };

    struct Energies
    {
        double ke_s = 0.0, ke_p = 0.0, ke_f = 0.0, u_b = 0.0;
    };

    struct Row
    {
        std::size_t step = 0;
        double t = 0.0, dt = 0.0, z_c = 0.0, v_z = 0.0;
        double KE_s = 0.0, KE_p = 0.0, KE_f = 0.0, U_b = 0.0, U_int = 0.0, E_tot = 0.0;
        double Fs[3] = {}, Ff[3] = {};
        std::size_t n_bonds = 0, n_broken_v = 0, n_broken_h = 0;
        double w_f = 0.0, w_p = 0.0;
        std::vector<std::size_t> pv, ph; // по слоям (только при n_plies >= 2)
        std::vector<double> wp;          // прогиб по слоям (только при n_plies >= 2)
    };

    struct ContactFF
    {
        std::size_t count = 0;
        double delta_d = 0.0;
        // Только при n_plies >= 2: «внутри слоя» и «между слоями»; max проникания / d.
        std::size_t intra = 0, inter = 0;
        double delta_intra = 0.0, delta_inter = 0.0;
    };

    // Свёртка по примитивам контакта: целые счётчики и max (детерминирована).
    struct SplitAcc
    {
        std::size_t intra = 0, inter = 0;
        double dmax_intra = 0.0, dmax_inter = 0.0;
    };

    void add_clamped(bool fiber, std::size_t i, std::uint8_t mask)
    {
        Clamped c{fiber, i, mask, {0, 0, 0}};
        for (int d = 0; d < 3; ++d) c.x0[d] = fiber ? m_.fiber.x[d][i] : m_.solid.x[d][i];
        clamped_.push_back(c);
    }

    double sum_mass_solid() const
    {
        double m = 0.0;
        for (std::size_t i = 0; i < nt_; ++i) m += m_.solid.mass[i];
        return m;
    }
    double sum_mass_fiber() const
    {
        double m = 0.0;
        for (std::size_t k = 0; k < nb_; ++k) m += m_.fiber.mass[k];
        return m;
    }

    double ke_range(std::size_t i0, std::size_t i1) const
    {
        const ParticleSystem& s = m_.solid;
        double e = 0.0;
        for (std::size_t i = i0; i < i1; ++i)
            for (int d = 0; d < 3; ++d) e += 0.5 * s.mass[i] * s.vx[d][i] * s.vx[d][i];
        return e;
    }

    Energies energies() const
    {
        Energies e;
        e.ke_s = ke_range(sc_.ball_begin, sc_.ball_end);
        e.ke_p = ke_range(sc_.plate_begin, sc_.plate_end);
        e.ke_f = clamp_e2e::kinetic_fiber(m_.fiber);
        // bond_energy: sum 1/2 (EA/l0)(l - l0)^2 по растянутым связям = sum 1/2 E A l0 eps^2
        e.u_b = clamp_e2e::bond_energy(m_);
        return e;
    }

    void count_broken(std::size_t& nv, std::size_t& nh)
    {
        std::fill(cnt_.begin(), cnt_.end(), 0u);
        for (const auto& b : m_.topology.bonds)
        {
            const int y = yarn_of_[b.i];
            if (y >= 0) ++cnt_[static_cast<std::size_t>(y)];
        }
        nv = nh = 0;
        if (nply_ >= 2)
        {
            std::fill(brk_pv_.begin(), brk_pv_.end(), std::size_t{0});
            std::fill(brk_ph_.begin(), brk_ph_.end(), std::size_t{0});
        }
        for (std::size_t y = 0; y < cnt_.size(); ++y)
            if (cnt_[y] < orig_bonds_[y])
            {
                const bool v = sc_.yarns[y].direction == 'y';
                if (v) ++nv;
                else ++nh;
                if (nply_ >= 2) ++(v ? brk_pv_ : brk_ph_)[static_cast<std::size_t>(sc_.yarns[y].ply - 1)];
            }
    }

    // Fiber-Fiber: принятые контакты и max проникание/d (боевой поиск и диагностика ядра).
    ContactFF contact_ff() const
    {
        ContactFF r;
        if (nb_ == 0) return r;
        // Постоянные члены монитора вместо локальных объектов на каждый вызов.
        interactions::build_fiber_contact_search(diag_executor_, m_.fiber, m_.topology, m_.materials, diag_search_);
        if (nply_ >= 2) return split_contacts();
        const auto dg = interactions::compute_fiber_contact_diagnostics(diag_executor_, m_.fiber, m_.materials, diag_search_);
        r.count = dg.contact_count;
        r.delta_d = d_fiber_ > 0.0 ? dg.max_penetration / d_fiber_ : 0.0;
        return r;
    }

    // Послойный обход diag_search_.primitives (по образцу compute_fiber_contact_diagnostics,
    // include/interactions/fiber_contact.hpp): принятый контакт -- внутри слоя, если ply_of_[a0] == ply_of_[b0].
    ContactFF split_contacts() const
    {
        ContactFF r;
        const std::size_t n = diag_search_.primitives.size();
        const fiber::FiberState& fb = m_.fiber;
        const std::vector<int>& po = ply_of_;
        const SplitAcc acc = diag_executor_.parallel_reduce(
            n, SplitAcc{},
            [&](std::size_t k) -> SplitAcc {
                SplitAcc a;
                const interactions::FiberContactPrimitive& prim = diag_search_.primitives[k];
                const interactions::FiberContactGeometry geom =
                    interactions::compute_fiber_contact_geometry(fb, m_.materials, prim);
                if (geom.rejected_by_ownership || geom.degenerate || !geom.contact) return a;
                if (po[prim.a0] == po[prim.b0])
                {
                    a.intra = 1;
                    a.dmax_intra = geom.delta;
                }
                else
                {
                    a.inter = 1;
                    a.dmax_inter = geom.delta;
                }
                return a;
            },
            [](SplitAcc x, SplitAcc y) {
                SplitAcc z;
                z.intra = x.intra + y.intra;
                z.inter = x.inter + y.inter;
                z.dmax_intra = std::max(x.dmax_intra, y.dmax_intra);
                z.dmax_inter = std::max(x.dmax_inter, y.dmax_inter);
                return z;
            });
        r.intra = acc.intra;
        r.inter = acc.inter;
        r.count = acc.intra + acc.inter;
        r.delta_intra = d_fiber_ > 0.0 ? acc.dmax_intra / d_fiber_ : 0.0;
        r.delta_inter = d_fiber_ > 0.0 ? acc.dmax_inter / d_fiber_ : 0.0;
        r.delta_d = std::max(r.delta_intra, r.delta_inter);
        return r;
    }

    // ---- пакет слоёв (ТЗ 2T/4T §4.3) -------------------------------------
    void init_plies()
    {
        nply_ = sc_.has_fabric ? sc_.n_plies : 1;
        if (nply_ < 2) return;
        const std::size_t N = static_cast<std::size_t>(nply_);
        ply_of_.assign(nb_, 0);
        for (std::size_t k = 0; k < N; ++k)
            for (std::size_t b = sc_.ply_layout[k].bead_begin; b < sc_.ply_layout[k].bead_end && b < nb_; ++b)
                ply_of_[b] = static_cast<int>(k);
        // вершинные бусины пересечений (как в ply_vertex_gaps): у основы j вершина i, у утка i вершина j
        const int nx = sc_.W.params.nx, ny = sc_.W.params.ny;
        vb_.assign(N, {});
        for (std::size_t k = 0; k < N; ++k)
        {
            vb_[k].reserve(2 * static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
            for (int j = 0; j < nx; ++j)
                for (int i = 0; i < ny; ++i)
                {
                    const std::size_t lx = static_cast<std::size_t>(j), ly = static_cast<std::size_t>(nx + i);
                    const std::size_t yx = sc_.ply_layout[k].yarn_begin + lx, yy = sc_.ply_layout[k].yarn_begin + ly;
                    vb_[k].push_back(sc_.yarns[yx].bead_begin + sc_.plies[k].fibers[lx].vertex_bead[static_cast<std::size_t>(i)]);
                    vb_[k].push_back(sc_.yarns[yy].bead_begin + sc_.plies[k].fibers[ly].vertex_bead[static_cast<std::size_t>(j)]);
                }
        }
        // концы свободных нитей (вдоль x): стартовые x
        end_x0_.assign(sc_.yarns.size(), {0.0, 0.0});
        for (std::size_t y = 0; y < sc_.yarns.size(); ++y)
            if (sc_.yarns[y].direction == 'x')
            {
                end_x0_[y][0] = m_.fiber.x[0][sc_.yarns[y].bead_begin];
                end_x0_[y][1] = m_.fiber.x[0][sc_.yarns[y].bead_end - 1];
            }
        brk_pv_.assign(N, 0);
        brk_ph_.assign(N, 0);
        pullin_cur_.assign(N, 0.0);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        st_.n_plies = nply_;
        st_.n_broken_v_p.assign(N, 0);
        st_.n_broken_h_p.assign(N, 0);
        st_.t_first_break_p.assign(N, -1.0);
        st_.w_f_at_first_break_p.assign(N, nan);
        st_.pullin_at_first_break_p.assign(N, nan);
        st_.pullin_end_p.assign(N, 0.0);
        st_.w_f_end_p.assign(N, 0.0);
        const PlyGapReport pg = ply_vertex_gaps(sc_); // по начальной укладке: монитор создаётся до прогона
        st_.ply_gap_min = pg.min_gap;
        st_.ply_gap_max = pg.max_gap;
    }

    // Послойные величины шага: прогиб по слоям, первый разрыв слоя, втягивание кромок.
    void plies_step(Row& r, double time)
    {
        const fiber::FiberState& fb = m_.fiber;
        const std::size_t N = static_cast<std::size_t>(nply_);
        r.wp.assign(N, 0.0);
        for (std::size_t k = 0; k < N; ++k)
        {
            double w = 0.0;
            for (std::size_t b = sc_.ply_layout[k].bead_begin; b < sc_.ply_layout[k].bead_end; ++b)
                w = std::max(w, z0_bead_[b] - fb.x[2][b]);
            r.wp[k] = w;
        }
        // втягивание: Delta = |x_0| - |x| концевых бусин свободных нитей слоя, максимум по нитям
        for (std::size_t k = 0; k < N; ++k)
        {
            double pm = -std::numeric_limits<double>::infinity();
            for (std::size_t y = sc_.ply_layout[k].yarn_begin; y < sc_.ply_layout[k].yarn_end; ++y)
            {
                if (sc_.yarns[y].direction != 'x') continue;
                pm = std::max(pm, std::fabs(end_x0_[y][0]) - std::fabs(fb.x[0][sc_.yarns[y].bead_begin]));
                pm = std::max(pm, std::fabs(end_x0_[y][1]) - std::fabs(fb.x[0][sc_.yarns[y].bead_end - 1]));
            }
            pullin_cur_[k] = std::isfinite(pm) ? pm : 0.0;
        }
        for (std::size_t k = 0; k < N; ++k)
            if (st_.t_first_break_p[k] < 0.0 && (r.pv[k] + r.ph[k]) > 0)
            {
                st_.t_first_break_p[k] = time;
                st_.w_f_at_first_break_p[k] = r.wp[k];
                st_.pullin_at_first_break_p[k] = pullin_cur_[k];
            }
    }

    // Нарушения порядка слоёв: пары (k, k+1) и пересечения, где средняя z слоя k <= средней z слоя k+1.
    std::size_t count_ply_order_violations() const
    {
        const fiber::FiberState& fb = m_.fiber;
        std::size_t nv = 0;
        for (std::size_t k = 0; k + 1 < vb_.size(); ++k)
        {
            const std::vector<std::size_t>& a = vb_[k];
            const std::vector<std::size_t>& b = vb_[k + 1];
            for (std::size_t q = 0; q + 1 < a.size(); q += 2)
            {
                const double za = 0.5 * (fb.x[2][a[q]] + fb.x[2][a[q + 1]]);
                const double zb = 0.5 * (fb.x[2][b[q]] + fb.x[2][b[q + 1]]);
                if (za <= zb) ++nv;
            }
        }
        return nv;
    }

    struct SolidGrid
    {
        double org[3] = {0, 0, 0};
        double W = 0.0;
        long long cell(double x, int d) const { return static_cast<long long>(std::floor((x - org[d]) / W)); }
    };

    // Дорогая диагностика: ff (боевое ядро), fs и туннелирование (сетка), m_ej.
    void diagnostics(bool first_step)
    {
        const ParticleSystem& s = m_.solid;
        const fiber::FiberState& fb = m_.fiber;

        // m_ej: частицы пластины дальше чем на t за исходной тыльной гранью и с v_z < 0
        if (sc_.has_plate)
        {
            const double zlim = sc_.plate_z_bot - sc_.plate_t;
            double mej = 0.0;
            for (std::size_t i = sc_.plate_begin; i < sc_.plate_end; ++i)
                if (s.x[2][i] < zlim && s.vx[2][i] < 0.0) mej += s.mass[i];
            st_.m_ej = mej;
            st_.m_ej_max = std::max(st_.m_ej_max, mej);
        }
        if (nb_ == 0) return;

        // Fiber-Fiber
        {
            const ContactFF c = contact_ff();
            st_.ff_contacts_max = std::max(st_.ff_contacts_max, c.count);
            st_.max_delta_ff_d = std::max(st_.max_delta_ff_d, c.delta_d);
            // При N >= 2 contacts_per_crossing -- контакты внутри слоя на N nx ny (при N = 1 то же, что раньше)
            const double cpc = crossings_ > 0.0
                                   ? static_cast<double>(nply_ >= 2 ? c.intra : c.count) / (crossings_ * (nply_ >= 2 ? nply_ : 1))
                                   : 0.0;
            st_.contacts_per_crossing = std::max(st_.contacts_per_crossing, cpc);
            double cpi = 0.0;
            if (nply_ >= 2)
            {
                cpi = crossings_ > 0.0 ? static_cast<double>(c.inter) / (crossings_ * (nply_ - 1)) : 0.0;
                st_.contacts_per_crossing_inter = std::max(st_.contacts_per_crossing_inter, cpi);
                st_.max_delta_ff_inter_d = std::max(st_.max_delta_ff_inter_d, c.delta_inter);
                st_.ff_intra_last = c.intra;
                st_.ff_inter_last = c.inter;
                const std::size_t pov = count_ply_order_violations();
                st_.ply_order_violations = std::max(st_.ply_order_violations, pov);
                if (pov > 0 && st_.t_ply_order_first < 0.0) st_.t_ply_order_first = st_.t_end;
                if (first_step) st_.ply_order_violations_step1 = pov;
            }
            if (first_step)
            {
                st_.ff_contacts_step1 = c.count;
                st_.ff_delta_d_step1 = c.delta_d;
                st_.contacts_per_crossing_step1 = cpc;
                st_.ff_intra_step1 = c.intra;
                st_.ff_inter_step1 = c.inter;
                st_.contacts_per_crossing_inter_step1 = cpi;
            }
        }

        // Fiber-Solid и туннелирование: сегменты нити против Solid-частиц, равномерная сетка
        const auto& bonds = m_.topology.bonds;
        if (bonds.empty() || nt_ == 0) return;
        double rs_max = 0.0;
        SolidGrid g;
        g.org[0] = g.org[1] = g.org[2] = 1e300;
        for (std::size_t k = 0; k < nt_; ++k)
        {
            rs_max = std::max(rs_max, interactions::solid_contact_radius(s, k));
            for (int d = 0; d < 3; ++d) g.org[d] = std::min(g.org[d], s.x[d][k]);
        }
        double rf_max = 0.0;
        for (const auto& mat : m_.materials) rf_max = std::max(rf_max, 0.5 * mat.d);
        g.W = rf_max + rs_max;
        if (!(g.W > 0.0)) return;
        // Ячейки сетки (ключ ячейки, индекс частицы), отсортированы; буфер постоянный.
        grid_cells_.clear();
        grid_cells_.reserve(nt_);
        for (std::size_t k = 0; k < nt_; ++k)
        {
            const std::uint64_t key = detail::cell_key(g.cell(s.x[0][k], 0), g.cell(s.x[1][k], 1), g.cell(s.x[2][k], 2));
            grid_cells_.emplace_back(key, static_cast<std::uint32_t>(k));
        }
        std::sort(grid_cells_.begin(), grid_cells_.end());

        // Метки туннелирования; буфер постоянный.
        tun_.clear();
        tun_.resize(nb_, 0);
        double dmax_fs = 0.0;
        for (const auto& b : bonds)
        {
            double a[3], bb[3], lo[3], hi[3];
            for (int d = 0; d < 3; ++d)
            {
                a[d] = fb.x[d][b.i];
                bb[d] = fb.x[d][b.j];
                lo[d] = std::min(a[d], bb[d]) - g.W;
                hi[d] = std::max(a[d], bb[d]) + g.W;
            }
            const double rf = 0.5 * m_.materials[static_cast<std::size_t>(fb.material_id[b.i])].d;
            const long long cx0 = g.cell(lo[0], 0), cx1 = g.cell(hi[0], 0);
            const long long cy0 = g.cell(lo[1], 1), cy1 = g.cell(hi[1], 1);
            const long long cz0 = g.cell(lo[2], 2), cz1 = g.cell(hi[2], 2);
            for (long long cx = cx0; cx <= cx1; ++cx)
                for (long long cy = cy0; cy <= cy1; ++cy)
                    for (long long cz = cz0; cz <= cz1; ++cz)
                    {
                        const std::uint64_t key = detail::cell_key(cx, cy, cz);
                        auto it = std::lower_bound(grid_cells_.begin(), grid_cells_.end(),
                                                   std::make_pair(key, static_cast<std::uint32_t>(0)));
                        for (; it != grid_cells_.end() && it->first == key; ++it)
                        {
                            const std::size_t k = it->second;
                            const double p[3] = {s.x[0][k], s.x[1][k], s.x[2][k]};
                            double q[3];
                            const double dist = detail::point_segment_dist(p, a, bb, q);
                            const double rs = interactions::solid_contact_radius(s, k);
                            const double delta = rf + rs - dist;
                            if (delta > 0.0) dmax_fs = std::max(dmax_fs, delta / (2.0 * rf));
                            // Глубокое проникание в шар: ближе r_f + r_s - 0.5 r_f (ТЗ Guo 2019 §6.3, n_tunnel)
                            if (k < sc_.ball_end && dist < rf + rs - 0.5 * rf)
                            {
                                tun_[b.i] = 1;
                                tun_[b.j] = 1;
                            }
                        }
                    }
        }
        st_.max_delta_fs_d = std::max(st_.max_delta_fs_d, dmax_fs);
        std::size_t nt = 0;
        for (char c : tun_) nt += (c != 0);
        st_.n_tunnel = std::max(st_.n_tunnel, nt);
    }

    // Критерий исхода (ТЗ Guo 2019 §6.2).
    void update_outcome(double t, double vz, double zc)
    {
        if (st_.outcome != Outcome::Running) return;
        hist_.emplace_back(t, vz);
        while (hist_.size() > 1 && hist_[1].first <= t - kOutcomeWindow) hist_.pop_front();
        if (vz >= 0.0)
        {
            if (t_pos_ < 0.0) t_pos_ = t;
        }
        else
            t_pos_ = -1.0;

        const bool full = hist_.front().first <= t - kOutcomeWindow;
        if (full && vz < 0.0 && std::fabs(vz - hist_.front().second) < kOutcomeDv && zc < sc_.z_b - sc_.R_eff)
        {
            st_.outcome = Outcome::Perforated;
            st_.V_R = -vz;
        }
        else if (t_pos_ >= 0.0 && t - t_pos_ >= kOutcomeWindow)
        {
            st_.outcome = Outcome::NotPerforated;
            st_.V_R = 0.0;
        }
        else if (t >= opt_.t_max)
        {
            st_.outcome = Outcome::Undetermined;
        }
        if (st_.outcome != Outcome::Running)
        {
            st_.t_outcome = t;
            if (opt_.stop_on_outcome && opt_.guard)
            {
                request_stop(*opt_.guard);
                st_.stopped_by_monitor = true;
            }
        }
    }

    void write_row(const Row& r)
    {
        if (!csv_) return;
        std::fprintf(csv_,
                     "%zu,%.12e,%.9e,%.12e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.12e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%zu,%zu,%zu,%.9e,%.9e%s",
                     r.step, r.t, r.dt, r.z_c, r.v_z, r.KE_s, r.KE_p, r.KE_f, r.U_b, r.U_int, r.E_tot, r.Fs[0], r.Fs[1],
                     r.Fs[2], r.Ff[0], r.Ff[1], r.Ff[2], r.n_bonds, r.n_broken_v, r.n_broken_h, r.w_f, r.w_p,
                     (nply_ >= 2 || proj_on_) ? "" : "\n");
        if (nply_ >= 2)
        {
            for (std::size_t k = 0; k < r.pv.size(); ++k) std::fprintf(csv_, ",%zu", r.pv[k]);
            for (std::size_t k = 0; k < r.ph.size(); ++k) std::fprintf(csv_, ",%zu", r.ph[k]);
            for (std::size_t k = 0; k < r.wp.size(); ++k) std::fprintf(csv_, ",%.9e", r.wp[k]);
            std::fprintf(csv_, "%s", proj_on_ ? "" : "\n");
        }
        if (proj_on_)
        {
            // Окно -- с прошлой записанной строки; пустое окно (нет шагов с проекцией) -- нули, warm_frac = 1.
            std::fprintf(csv_, ",%zu,%zu,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n", pw_.n_ff, pw_.n_fs, pw_.max_ck, pw_.max_cs,
                         pw_.resid, pw_.resid_rel, pw_.warm_frac, pw_.W_pos);
            pw_ = ProjWindow{};
        }
    }

    // Окно накопления столбцов проекции между записанными строками monitor.csv.
    //   n_ff, n_fs, max_ck, max_cs, resid, resid_rel -- максимумы по шагам окна;
    //   warm_frac -- МИНИМУМ по шагам окна с активными слотами (1, если слотов не было);
    //   W_pos     -- max(0, max W_stage по шагам окна): положительная работа стадии (теоретически 0).
    struct ProjWindow
    {
        std::size_t n_ff = 0, n_fs = 0;
        double max_ck = 0.0, max_cs = 0.0, resid = 0.0, resid_rel = 0.0;
        double warm_frac = 1.0;
        double W_pos = 0.0;
    };

    void accumulate_proj(const interactions::NormalProjectionStats& ps)
    {
        pw_.n_ff = std::max(pw_.n_ff, ps.n_slots_ff);
        pw_.n_fs = std::max(pw_.n_fs, ps.n_slots_fs);
        pw_.max_ck = std::max(pw_.max_ck, ps.max_ck);
        pw_.max_cs = std::max(pw_.max_cs, ps.max_cs);
        pw_.resid = std::max(pw_.resid, ps.resid);
        pw_.resid_rel = std::max(pw_.resid_rel, ps.resid_rel);
        const std::size_t n_slots = ps.n_slots_ff + ps.n_slots_fs;
        if (n_slots > 0) pw_.warm_frac = std::min(pw_.warm_frac, ps.warm_frac);
        pw_.W_pos = std::max(pw_.W_pos, ps.W_stage);
        // весь прогон
        st_.proj_resid_max = std::max(st_.proj_resid_max, ps.resid);
        st_.proj_resid_rel_max = std::max(st_.proj_resid_rel_max, ps.resid_rel);
        st_.W_proj_pos_max = std::max(st_.W_proj_pos_max, ps.W_stage);
        st_.proj_slots_max = std::max(st_.proj_slots_max, n_slots);
    }

    static constexpr std::uint64_t kStateMagic = 0x3130'4e4f'4d4f'5547ull; // "GUOMON01"

    // Все поля, меняющиеся после конструктора (накопленные величины, исход, окно проекции, последняя строка).
    // Величины конструктора (x0 зажатых, массы, z0, U_int0, пакет слоёв) -- от исходной сцены, она та же.
    // ДОБАВЛЯЯ поле в MonitorStats или изменяемый член класса -- добавить его сюда.
    template <class F>
    void state_io(F& f)
    {
        MonitorStats& s = st_;
        f(s.steps), f(s.t_end), f(s.dt_first), f(s.dt_min);
        f(s.E_tot_0), f(s.E_tot_max), f(s.E_tot_end);
        f(s.KE_s_0), f(s.KE_s_end), f(s.KE_p_end), f(s.KE_f_end), f(s.U_b_end), f(s.U_int_end);
        f(s.z_c_end), f(s.v_z_end);
        f(s.F_s_z_max), f(s.F_s_z_max_t), f(s.F_f_z_max), f(s.F_f_z_max_t), f(s.n_clamped_s), f(s.n_clamped_f);
        f(s.n_bonds_0), f(s.n_bonds_end), f(s.n_broken_v), f(s.n_broken_h);
        f(s.t_first_break), f(s.w_f_at_first_break), f(s.w_f_max), f(s.w_p_max), f(s.w_f_end);
        f(s.eps_p_max), f(s.W_p_total), f(s.n_bonds_plastic);
        f(s.n_plies), f(s.n_broken_v_p), f(s.n_broken_h_p), f(s.t_first_break_p), f(s.w_f_at_first_break_p);
        f(s.pullin_at_first_break_p), f(s.pullin_end_p), f(s.w_f_end_p);
        f(s.contacts_per_crossing_inter), f(s.contacts_per_crossing_inter_step1), f(s.max_delta_ff_inter_d);
        f(s.ff_intra_initial), f(s.ff_inter_initial), f(s.ff_intra_step1), f(s.ff_inter_step1);
        f(s.ff_intra_last), f(s.ff_inter_last), f(s.ply_order_violations), f(s.ply_order_violations_step1);
        f(s.t_ply_order_first), f(s.ply_gap_min), f(s.ply_gap_max);
        f(s.m_ej), f(s.m_ej_max), f(s.n_tunnel), f(s.max_delta_ff_d), f(s.max_delta_fs_d);
        f(s.contacts_per_crossing), f(s.contacts_per_crossing_step1), f(s.ff_contacts_initial), f(s.ff_delta_d_initial);
        f(s.ff_contacts_step1), f(s.ff_delta_d_step1), f(s.ff_contacts_max);
        f(s.plate_vmax_step1), f(s.bead_vmax_step1);
        f(s.clamp_violations), f(s.clamp_first_violation), f(s.mass_violations);
        f(s.outcome), f(s.V_R), f(s.t_outcome), f(s.stopped_by_monitor);
        f(s.monitor_wall_s);
        f(s.proj_resid_max), f(s.proj_resid_rel_max), f(s.W_proj_pos_max), f(s.proj_slots_max);

        f(last_bonds_), f(broken_v_), f(broken_h_), f(cnt_), f(brk_pv_), f(brk_ph_), f(pullin_cur_);
        f(t_pos_), f(pw_), f(last_written_);
        Row& r = last_;
        f(r.step), f(r.t), f(r.dt), f(r.z_c), f(r.v_z), f(r.KE_s), f(r.KE_p), f(r.KE_f), f(r.U_b), f(r.U_int), f(r.E_tot);
        f(r.Fs), f(r.Ff), f(r.n_bonds), f(r.n_broken_v), f(r.n_broken_h), f(r.w_f), f(r.w_p), f(r.pv), f(r.ph), f(r.wp);

        // окно исхода (deque) и пластическая работа по связям (unordered_map, в файле -- по возрастанию ключа)
        std::vector<double> h;
        std::vector<std::uint64_t> keys;
        std::vector<double> vals;
        if constexpr (std::is_same_v<F, mon_detail::StateOut>)
        {
            for (const auto& e : hist_) h.push_back(e.first), h.push_back(e.second);
            for (const auto& kv : wp_bond_) keys.push_back(kv.first);
            std::sort(keys.begin(), keys.end());
            for (std::uint64_t k : keys) vals.push_back(wp_bond_.at(k));
        }
        f(h), f(keys), f(vals);
        if constexpr (!std::is_same_v<F, mon_detail::StateOut>)
        {
            if (h.size() % 2 != 0 || keys.size() != vals.size())
                throw std::runtime_error("состояние монитора в чекпоинте повреждено");
            hist_.clear();
            for (std::size_t i = 0; i < h.size(); i += 2) hist_.emplace_back(h[i], h[i + 1]);
            wp_bond_.clear();
            for (std::size_t i = 0; i < keys.size(); ++i) wp_bond_[keys[i]] = vals[i];
        }
    }

    // monitor.csv после рестарта с состоянием: заголовок и строки с шагом <= st_.steps (полные: столько же
    // полей, сколько в заголовке) сохраняются, дальше файл дописывается. Нет файла -- заводится заново.
    void reopen_csv_after_restart()
    {
        csv_pending_ = false;
        std::ifstream in(opt_.csv_path);
        if (!in)
        {
            std::fprintf(stderr, "ПРЕДУПРЕЖДЕНИЕ: монитор: нет %s -- строки до шага рестарта потеряны.\n",
                         opt_.csv_path.c_str());
            open_csv_fresh();
            return;
        }
        std::vector<std::string> keep;
        std::string line;
        std::size_t ncomma = 0, last_step = 0, dropped_tail = 0;
        while (std::getline(in, line))
        {
            if (keep.empty())
            {
                ncomma = static_cast<std::size_t>(std::count(line.begin(), line.end(), ','));
                keep.push_back(line);
                continue;
            }
            char* end = nullptr;
            const unsigned long long st = std::strtoull(line.c_str(), &end, 10);
            const bool full = end != line.c_str() && *end == ',' &&
                              static_cast<std::size_t>(std::count(line.begin(), line.end(), ',')) == ncomma;
            if (full && st <= st_.steps)
            {
                keep.push_back(line);
                last_step = static_cast<std::size_t>(st);
            }
            else
                ++dropped_tail;
        }
        in.close();
        const std::string tmp = opt_.csv_path + ".tmp";
        {
            std::ofstream out(tmp, std::ios::trunc);
            for (const auto& l : keep) out << l << '\n';
        }
        std::rename(tmp.c_str(), opt_.csv_path.c_str());
        csv_ = std::fopen(opt_.csv_path.c_str(), "a");
        std::fprintf(stderr, "монитор: состояние восстановлено на шаге %zu; monitor.csv: оставлено строк %zu (последняя -- шаг %zu), "
                             "отброшено после чекпоинта %zu.\n", st_.steps, keep.empty() ? 0 : keep.size() - 1, last_step,
                     dropped_tail);
    }

    const GuoScene& sc_;
    const fiber::Model& m_;
    MonitorOptions opt_;
    std::size_t nt_, nb_;
    MonitorStats st_;

    std::vector<Clamped> clamped_;
    double mass_s0_ = 0.0, mass_f0_ = 0.0, ball_mass_ = 0.0, U_int0_ = 0.0;
    std::vector<double> z0_bead_, z0_back_;
    std::vector<int> yarn_of_;
    std::vector<std::size_t> orig_bonds_;
    std::vector<std::size_t> cnt_;
    std::size_t last_bonds_ = 0, broken_v_ = 0, broken_h_ = 0;
    double d_fiber_ = 0.0, crossings_ = 0.0;

    // Пакет слоёв (n_plies >= 2)
    int nply_ = 1;
    std::vector<int> ply_of_;                       // бусина -> слой (с нуля)
    std::vector<std::vector<std::size_t>> vb_;      // по слоям: пары вершинных бусин пересечений (основа, уток)
    std::vector<std::array<double, 2>> end_x0_;     // по нитям: стартовые x концевых бусин
    std::vector<std::size_t> brk_pv_, brk_ph_;      // порванные нити по слоям
    std::vector<double> pullin_cur_;                // втягивание кромок по слоям на текущем шаге

    std::deque<std::pair<double, double>> hist_;
    double t_pos_ = -1.0;

    std::FILE* csv_ = nullptr;
    bool csv_pending_ = false; // restart: monitor.csv ещё не открыт (ждём restore_state или первого шага)
    bool restored_ = false;    // состояние восстановлено из чекпоинта
    bool proj_on_ = false; // режим проекции нормального контакта включён в конфиге сцены
    ProjWindow pw_;
    Row last_;
    bool last_written_ = true;

    // Персистентные буферы диагностики: переиспользуются в diagnostics() вместо
    // повторной аллокации (контакты FF, FS, туннелирование, m_ej).
    mutable rhs::CpuExecutor diag_executor_{1};
    mutable interactions::FiberContactSearchResult diag_search_;
    mutable std::vector<std::pair<std::uint64_t, std::uint32_t>> grid_cells_;
    mutable std::vector<char> tun_;
};

} // namespace guo2019
