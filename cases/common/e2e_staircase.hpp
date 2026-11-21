#pragma once
// e2e_staircase.hpp
// Выделение «лесенки» на кривой монитора: ступени (плато), провалы между
// ними, отношения соседних ступеней, период и возврат после провала.
//
// Новый код: у BC2 профиль — одно плато, и такой метрики не было. Кейсу BC6
// она нужна, потому что из-за большой разницы импедансов измеренный профиль
// там — убывающая лесенка из трёх ступеней с двумя провалами, и §6
// постановки требует сверять именно её.
//
// ПОЧЕМУ ПРОВАЛ ОПРЕДЕЛЯЕТСЯ КАК ГРАНИЦА МЕЖДУ СТУПЕНЯМИ.
// Наивный поиск минимума в окне ±4 точки на приборном эталоне BC6 даёт
// ДВЕНАДЦАТЬ кандидатов, из которых настоящих провала два: 114.93 м/с при
// 1.990 мкс и 105.79 м/с при 2.860 мкс.
//
// Первый отсев — по ВЫРАЖЕННОСТИ: насколько минимум ниже максимума сглаженной
// кривой СЛЕВА в пределах одного номинального периода (nominal_period), плюс
// слабое условие отскока СПРАВА (min_rebound, тоже в пределах периода).
// ИСТОРИЯ (2026-11-08): прежде бралась min(макс. слева, макс. справа ДО КОНЦА
// ЗАПИСИ). Хвост расчёта BC6 после ~2 мкс невоспроизводим, и при низком хвосте
// провал 1 (≈1.53 мкс) отбрасывался; без провала ступень 1 считалась медианой
// всего хвоста (127.9 м/с вместо полки 182.8) — ложный провал ворот плато.
// Отсев убирает шум оцифровки
// внутри плато (глубина 1-3 м/с), но НЕ убирает главного самозванца: внутри
// плато 1 эталона есть просадка до 160.8 м/с при 1.390 мкс с выраженностью
// 10.75 м/с — БОЛЬШЕЙ, чем у настоящего второго провала (≈10.3 м/с). Отбор
// «два самых выраженных» на этом и ломается: он выбирает самозванца и
// отбрасывает настоящий провал.
//
// Поэтому решающий критерий другой и отвечает сути лесенки: провал — это
// ГРАНИЦА МЕЖДУ СТУПЕНЯМИ, и уровень после него обязан упасть. Сравниваются
// медианы в окнах до и после минимума:
//   просадка внутри плато 1 : 172.4 → 171.4, падение 0.6 %  — отвергнута;
//   провал 1                : 171.5 → 141.0, падение 17.8 % — принят;
//   провал 2                : 141.0 → 115.5, падение 18.1 % — принят.
// Разделение на два порядка, порог 8 % лежит посередине с огромным запасом.
// Он же согласован с самой постановкой: §6 ждёт отношения ступеней
// 0.82 ± 0.02, то есть падения на 18 ± 2 %.
//
// ПОЧЕМУ МЕДИАНА, А НЕ СРЕДНЕЕ. Уровень ступени берётся медианой: эталон
// оцифрован с шумом, а осевой монитор расчёта может нести приосевой шум
// (§4 постановки), и одиночный выброс не должен сдвигать уровень.
//
// ПОЧЕМУ ОКНА ПЛАТО СТРОЯТСЯ ОТ НАЙДЕННЫХ СОБЫТИЙ, А НЕ ОТ АБСОЛЮТНЫХ ВРЕМЁН
// §5.2. Расчёт вправе дать другой период лесенки — это одна из проверяемых
// величин. Окна, привязанные к абсолютной шкале эталона, при другом периоде
// поехали бы на провалы, и уровень ступени вышел бы заниженным по причине,
// не имеющей отношения к уровню.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "e2e_curve.hpp"
#include "e2e_shape.hpp"

namespace e2e
{

struct StaircaseParams
{
    // Откуда начинать поиск провалов: t90 переднего фронта плюс отступ.
    // Отступ нужен, чтобы выброс на фронте и его спад не выглядели провалом.
    double search_lag     = 5.0e-8;   // с
    // Окно скользящего среднего ПО ВРЕМЕНИ. Именно по времени, а не по числу
    // точек: шаг неравномерен и у эталона (оцифровка), и у расчёта (шаг CFL).
    double smooth_window  = 4.0e-8;   // с
    // Минимальная выраженность провала в долях пика. Грубый предварительный
    // отсев шума оцифровки; решающий критерий — падение уровня ниже.
    double min_prominence = 0.04;
    // Минимальный промежуток между двумя засчитанными провалами.
    double min_separation = 3.0e-7;   // с
    // Окна слева и справа от минимума, по которым проверяется падение уровня.
    double step_window    = 2.0e-7;   // с
    // Зазор между минимумом и окнами: сам провал в них попадать не должен.
    double step_guard     = 4.0e-8;   // с
    // Минимальное падение уровня через провал, в долях уровня ДО него.
    double min_step_drop  = 0.08;
    // Доля интервала, отрезаемая с каждого края при усреднении ступени.
    double plateau_trim   = 0.15;
    // Номинальный период лесенки. Нужен в двух местах: (а) полуширина ЛОКАЛЬНОГО
    // окна, в котором мерится выраженность провала (максимумы слева и справа
    // берутся не дальше, чем на период от минимума, и не за пределами области
    // поиска/записи); (б) окно ступени 1 от фронта, когда провалов не найдено.
    // Значение — BC6 (§6 постановки, bc6::kPeriod = 0.88 мкс); для кривых с
    // иным периодом передавать своё.
    //
    // ПОЧЕМУ ЛОКАЛЬНОЕ ОКНО. Выраженность «до конца записи» зависит от хвоста:
    // в расчёте хвост после ~2 мкс невоспроизводим от прогона к прогону.
    // Настоящий провал на границе ступеней 1 и 2 (BC6, ≈1.53 мкс) неглубок
    // относительно СЛЕДУЮЩЕЙ полки (≈142 против ≈146 м/с), и его выраженность
    // «min(слева, справа до конца записи)» держалась на далёком отскоке после
    // 2.5 мкс; опустился хвост — провал пропал, а плато 1 уходило в медиану
    // всего хвоста. Поэтому глубина мерится от ЛОКАЛЬНОГО максимума слева (не
    // дальше периода), а справа требуется лишь заметный отскок (min_rebound).
    double nominal_period = 0.88e-6;  // с
    // Минимальный отскок сглаженной кривой справа от минимума (в пределах
    // nominal_period), в долях пика. Отличает дно от края ступени: на чистой
    // ступеньке без провала справа от точки излома подъёма нет. Порог около
    // уровня шума оцифровки (1-3 м/с ≈ 1 % пика): решающий критерий всё равно
    // падение уровня через минимум.
    double min_rebound    = 0.01;
};

struct StaircaseFeatures
{
    bool        valid = false;   // передний фронт прослежен
    std::string reason;
    double peak    = 0.0;
    double t_peak  = 0.0;
    double t_front = 0.0;        // t10 переднего фронта
    double t_mid   = 0.0;        // t50 переднего фронта: на размазанном
                                 // фронте именно он сопоставим с
                                 // теоретическим моментом прихода

    size_t dips = 0;                       // 0, 1 или 2
    std::array<double, 2> t_dip{};
    std::array<double, 2> v_dip{};
    std::array<double, 2> prominence{};

    size_t plateaus = 0;                   // 1..3
    std::array<double, 3> plateau{};       // медиана уровня
    std::array<double, 3> plateau_t0{};
    std::array<double, 3> plateau_t1{};
    std::array<size_t, 3> plateau_n{};     // сколько точек попало в окно

    std::array<double, 2> ratio{};         // plateau[k+1] / plateau[k]
    std::array<double, 2> period{};        // t_dip0 - t_front, t_dip1 - t_dip0
    std::array<double, 2> pullback{};      // plateau[k+1] - v_dip[k]
};

namespace detail
{

// Скользящее среднее по временно́му окну. Возвращает вектор той же длины.
inline std::vector<double> smooth_by_time(const Curve& c, double window)
{
    std::vector<double> out(c.v.size(), 0.0);
    if (!(window > 0.0)) { out = c.v; return out; }
    const double half = 0.5 * window;
    size_t lo = 0, hi = 0;
    double sum = 0.0;
    for (size_t i = 0; i < c.t.size(); ++i)
    {
        while (hi < c.t.size() && c.t[hi] <= c.t[i] + half) { sum += c.v[hi]; ++hi; }
        while (lo < hi && c.t[lo] < c.t[i] - half)          { sum -= c.v[lo]; ++lo; }
        const size_t n = hi - lo;
        out[i] = (n > 0) ? sum / static_cast<double>(n) : c.v[i];
    }
    return out;
}

inline double median_in_window(const Curve& c, double t0, double t1, size_t& n_out)
{
    std::vector<double> vals;
    for (size_t i = 0; i < c.t.size(); ++i)
    {
        if (c.t[i] >= t0 && c.t[i] <= t1) vals.push_back(c.v[i]);
    }
    n_out = vals.size();
    if (vals.empty()) return 0.0;
    std::sort(vals.begin(), vals.end());
    const size_t m = vals.size() / 2;
    return (vals.size() % 2 == 1) ? vals[m] : 0.5 * (vals[m - 1] + vals[m]);
}

struct DipCandidate
{
    size_t index;
    double t;
    double v;
    double prominence;
    double drop;     // относительное падение уровня через провал
};

} // namespace detail

inline StaircaseFeatures extract_staircase(const Curve& c, const StaircaseParams& par = {})
{
    StaircaseFeatures f;

    const CurveFeatures cf = extract_features(c);
    if (!cf.valid)
    {
        f.reason = cf.reason;
        return f;
    }
    f.peak    = cf.peak;
    f.t_peak  = cf.t_peak;
    f.t_front = cf.t10;
    f.t_mid   = cf.t50;
    f.valid   = true;

    // Область поиска провалов.
    const double t_search = cf.t90 + par.search_lag;
    const double t_end    = c.t.back();
    size_t i0 = 0;
    while (i0 < c.t.size() && c.t[i0] < t_search) ++i0;
    if (i0 + 3 >= c.t.size())
    {
        // Записи после фронта почти нет: ступеней выделить не из чего, но
        // сам фронт прослежен — это не повод объявлять кривую негодной.
        f.plateaus = 0;
        return f;
    }

    const std::vector<double> sm = detail::smooth_by_time(c, par.smooth_window);

    const size_t n = c.t.size();

    std::vector<detail::DipCandidate> cand;
    for (size_t i = i0 + 1; i + 1 < n; ++i)
    {
        if (!(sm[i] <= sm[i - 1] && sm[i] <= sm[i + 1])) continue;

        // Глубина — от ЛОКАЛЬНОГО максимума слева (не дальше периода лесенки,
        // не раньше начала поиска); справа нужен лишь отскок: минимум обязан
        // быть дном, а не краем ступени. Хвост записи на отбор не влияет.
        double max_l = sm[i], max_r = sm[i];
        for (size_t j = i; j-- > i0;)
        {
            if (c.t[i] - c.t[j] > par.nominal_period) break;
            max_l = std::max(max_l, sm[j]);
        }
        for (size_t j = i + 1; j < n; ++j)
        {
            if (c.t[j] - c.t[i] > par.nominal_period) break;
            max_r = std::max(max_r, sm[j]);
        }
        const double prom = max_l - sm[i];
        if (!(prom >= par.min_prominence * f.peak)) continue;
        if (!(max_r - sm[i] >= par.min_rebound * f.peak)) continue;

        // Решающий критерий: уровень ПОСЛЕ минимума обязан быть заметно ниже
        // уровня ДО него, иначе это просадка внутри ступени, а не граница
        // между ступенями.
        size_t n_before = 0, n_after = 0;
        const double lvl_before = detail::median_in_window(
            c, c.t[i] - par.step_guard - par.step_window, c.t[i] - par.step_guard, n_before);
        const double lvl_after = detail::median_in_window(
            c, c.t[i] + par.step_guard, c.t[i] + par.step_guard + par.step_window, n_after);
        // Меньше трёх точек в окне — подтвердить ступень нечем. Так
        // отбрасываются минимумы у самого края записи, где «падения» ещё или
        // уже не видно.
        if (n_before < 3 || n_after < 3) continue;
        if (!(lvl_before > 0.0)) continue;
        const double drop = (lvl_before - lvl_after) / lvl_before;
        if (!(drop >= par.min_step_drop)) continue;

        cand.push_back({i, c.t[i], sm[i], prom, drop});
    }

    // Отбор: сначала по величине падения (именно оно определяет ступень),
    // затем разведение по времени.
    std::sort(cand.begin(), cand.end(),
              [](const detail::DipCandidate& a, const detail::DipCandidate& b)
              { return a.drop > b.drop; });
    std::vector<detail::DipCandidate> picked;
    for (const auto& k : cand)
    {
        bool too_close = false;
        for (const auto& p : picked)
        {
            if (std::fabs(k.t - p.t) < par.min_separation) { too_close = true; break; }
        }
        if (too_close) continue;
        picked.push_back(k);
        if (picked.size() == 2) break;
    }
    std::sort(picked.begin(), picked.end(),
              [](const detail::DipCandidate& a, const detail::DipCandidate& b)
              { return a.t < b.t; });

    // Уточнение по НЕсглаженной кривой: сглаживание нужно было только для
    // устойчивого поиска, а сообщать полагается измеренные значения.
    f.dips = picked.size();
    for (size_t k = 0; k < picked.size(); ++k)
    {
        double best_t = picked[k].t, best_v = c.v[picked[k].index];
        for (size_t i = 0; i < n; ++i)
        {
            if (std::fabs(c.t[i] - picked[k].t) > par.smooth_window) continue;
            if (c.v[i] < best_v) { best_v = c.v[i]; best_t = c.t[i]; }
        }
        f.t_dip[k] = best_t;
        f.v_dip[k] = best_v;
        f.prominence[k] = picked[k].prominence;
    }

    // Окна ступеней строятся от найденных событий.
    const double trim = par.plateau_trim;
    auto set_plateau = [&](size_t k, double a, double b, size_t min_n = 1) -> bool
    {
        if (!(b > a)) return false;
        size_t cnt = 0;
        const double med = detail::median_in_window(c, a, b, cnt);
        if (cnt < min_n) return false;
        f.plateau_t0[k] = a;
        f.plateau_t1[k] = b;
        f.plateau[k] = med;
        f.plateau_n[k] = cnt;
        f.plateaus = k + 1;
        return true;
    };

    // Ступень 1 от фронта: [t10 + trim·T, t10 + (1 − trim)·T], обрезанное
    // концом записи. Не зависит от хвоста, поэтому им не размывается: медиана
    // «всего хвоста» (прежнее поведение при dips == 0) смешивала ступени 1…3
    // и уходила из ворот 170…190 по причине, не связанной с полкой.
    auto plateau1_from_front = [&]()
    {
        const double T = par.nominal_period;
        set_plateau(0, f.t_front + trim * T, std::min(f.t_front + (1.0 - trim) * T, t_end));
    };

    if (f.dips == 0)
    {
        // Ступени 2 и 3 без найденных провалов не выдумываются.
        plateau1_from_front();
    }
    else
    {
        const double s0 = f.t_dip[0] - f.t_front;
        // Вырожденное окно (пустое или меньше трёх точек) — запасной вариант
        // от фронта.
        if (!set_plateau(0, f.t_front + trim * s0, f.t_dip[0] - trim * s0, 3))
        {
            plateau1_from_front();
        }
        if (f.dips == 1)
        {
            const double s1 = t_end - f.t_dip[0];
            set_plateau(1, f.t_dip[0] + trim * s1, t_end);
        }
        else
        {
            const double s1 = f.t_dip[1] - f.t_dip[0];
            set_plateau(1, f.t_dip[0] + trim * s1, f.t_dip[1] - trim * s1);
            const double s2 = t_end - f.t_dip[1];
            set_plateau(2, f.t_dip[1] + trim * s2, t_end);
        }
    }

    for (size_t k = 0; k + 1 < f.plateaus && k < 2; ++k)
    {
        if (f.plateau[k] > 0.0) f.ratio[k] = f.plateau[k + 1] / f.plateau[k];
    }
    if (f.dips >= 1) f.period[0] = f.t_dip[0] - f.t_front;
    if (f.dips >= 2) f.period[1] = f.t_dip[1] - f.t_dip[0];
    for (size_t k = 0; k < f.dips && k + 1 < f.plateaus; ++k)
    {
        f.pullback[k] = f.plateau[k + 1] - f.v_dip[k];
    }

    return f;
}

inline void print_staircase(const char* who, const StaircaseFeatures& f)
{
    constexpr double us = 1e6;
    if (!f.valid)
    {
        std::printf("  info : лесенка %s: не выделена (%s)\n", who, f.reason.c_str());
        return;
    }
    std::printf("  info : лесенка %s: пик %.2f м/с при %.4f мкс; фронт (10 %%) при %.4f мкс\n",
                who, f.peak, f.t_peak * us, f.t_front * us);
    for (size_t k = 0; k < f.plateaus; ++k)
    {
        std::printf("         ступень %zu: %.2f м/с (медиана по %zu точкам, окно %.4f…%.4f мкс)\n",
                    k + 1, f.plateau[k], f.plateau_n[k],
                    f.plateau_t0[k] * us, f.plateau_t1[k] * us);
    }
    for (size_t k = 0; k < f.dips; ++k)
    {
        std::printf("         провал %zu: %.2f м/с при %.4f мкс (выраженность %.2f м/с)\n",
                    k + 1, f.v_dip[k], f.t_dip[k] * us, f.prominence[k]);
    }
    if (f.plateaus >= 2) std::printf("         отношение ступеней 1→2: %.3f\n", f.ratio[0]);
    if (f.plateaus >= 3) std::printf("         отношение ступеней 2→3: %.3f\n", f.ratio[1]);
    if (f.dips >= 1) std::printf("         период фронт→провал 1: %.3f мкс\n", f.period[0] * us);
    if (f.dips >= 2) std::printf("         период провал 1→провал 2: %.3f мкс\n", f.period[1] * us);
    for (size_t k = 0; k < f.dips && k + 1 < f.plateaus; ++k)
    {
        std::printf("         возврат после провала %zu: %+.2f м/с\n", k + 1, f.pullback[k]);
    }
}

} // namespace e2e
