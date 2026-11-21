#pragma once
// e2e_shape.hpp
// Сравнение кривых ПО ФОРМЕ: признаки (пик, фронт 10/50/90 %, спад до 70 %),
// совмещение по t50 переднего фронта, безразмерные отношения и ворота.
// Плюс самопроверка самого компаратора.
//
// Происхождение — tests/test_e2e_bc2.cpp, раздел 2 и раздел 7b.
// Сам BC2 на этот заголовок не переключён, см. пояснение в e2e_harness.hpp.
//
// Что параметризовано против оригинала:
//   * пороги ворот — были четыре constexpr на уровне файла, стали ShapeGates;
//   * текст оговорки про сдвиг времени — был цитатой из Grady & Moody,
//     привязанной к BC2, стал аргументом print_shape_report();
//   * эталон для самопроверки — BC2 читал его сам, здесь он передаётся.
//
// ПОЧЕМУ СРАВНИВАЕТСЯ ФОРМА, А НЕ ЗНАЧЕНИЯ В ОБЩЕМ АБСОЛЮТНОМ ВРЕМЕНИ.
// В обоих кейсах (BC2 и BC6) начало отсчёта времени у приборного эталона
// произвольно: времена прихода первой волны в этих опытах не измерялись.
// Сравнение по абсолютному времени измеряло бы величину, лишённую смысла,
// поэтому кривые совмещаются по t50 переднего фронта, а сдвиг печатается
// справочно и в вердикт не входит.

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

#include "e2e_curve.hpp"
#include "e2e_harness.hpp"

namespace e2e
{

// Уровни в долях пика, по которым берутся признаки формы.
inline constexpr double kLevelLo   = 0.10;   // низ переднего фронта
inline constexpr double kLevelMid  = 0.50;   // опорная точка совмещения
inline constexpr double kLevelHi   = 0.90;   // верх переднего фронта
inline constexpr double kLevelFall = 0.70;   // уровень, до которого меряется спад

struct CurveFeatures
{
    bool        valid = false;
    std::string reason;
    double peak   = 0.0;   // максимум v, м/с
    double t_peak = 0.0;   // момент пика, с
    double t10    = 0.0;   // пересечение 10 % пика на переднем фронте, с
    double t50    = 0.0;   // 50 % — опорная точка совмещения, с
    double t90    = 0.0;   // 90 %, с
    double rise   = 0.0;   // длительность фронта 10-90 %, с
    double fall   = 0.0;   // от пика до спада на 70 % пика, с
    bool   fall_valid = false;   // спад до 70 % может не попасть в окно записи
};

// Момент, когда кривая, идя вверх, ПОСЛЕДНИЙ раз перед пиком пересекает level.
// Поиск от пика назад: на реальном сигнале уровень пересекается не один раз
// (шум, упругий предвестник), и содержательно последнее пересечение перед пиком.
inline bool front_crossing(const Curve& c, size_t i_peak, double level, double& t_out)
{
    for (size_t i = i_peak; i-- > 0;)
    {
        if (c.v[i] <= level && c.v[i + 1] >= level)
        {
            const double den = c.v[i + 1] - c.v[i];
            t_out = (den > 0.0)
                        ? c.t[i] + ((level - c.v[i]) / den) * (c.t[i + 1] - c.t[i])
                        : c.t[i + 1];
            return true;
        }
    }
    return false;
}

// Первое пересечение уровня вниз ПОСЛЕ пика.
inline bool fall_crossing(const Curve& c, size_t i_peak, double level, double& t_out)
{
    for (size_t i = i_peak; i + 1 < c.v.size(); ++i)
    {
        if (c.v[i] >= level && c.v[i + 1] <= level)
        {
            const double den = c.v[i] - c.v[i + 1];
            t_out = (den > 0.0)
                        ? c.t[i] + ((c.v[i] - level) / den) * (c.t[i + 1] - c.t[i])
                        : c.t[i + 1];
            return true;
        }
    }
    return false;
}

inline CurveFeatures extract_features(const Curve& c)
{
    CurveFeatures f;
    if (c.t.size() < 2)
    {
        f.reason = "кривая содержит меньше двух точек";
        return f;
    }
    for (size_t i = 0; i < c.t.size(); ++i)
    {
        if (!std::isfinite(c.t[i]) || !std::isfinite(c.v[i]))
        {
            std::ostringstream os;
            os << "в кривой есть неконечное значение (точка #" << i << ")";
            f.reason = os.str();
            return f;
        }
    }
    size_t i_peak = 0;
    for (size_t i = 1; i < c.v.size(); ++i)
    {
        if (c.v[i] > c.v[i_peak]) i_peak = i;
    }
    f.peak   = c.v[i_peak];
    f.t_peak = c.t[i_peak];
    if (!(f.peak > 0.0))
    {
        f.reason = "у кривой нет положительного пика (сигнал не пришёл?)";
        return f;
    }
    if (i_peak == 0)
    {
        f.reason = "пик пришёлся на первую точку — переднего фронта в записи нет";
        return f;
    }
    if (!front_crossing(c, i_peak, kLevelLo * f.peak, f.t10) ||
        !front_crossing(c, i_peak, kLevelMid * f.peak, f.t50) ||
        !front_crossing(c, i_peak, kLevelHi * f.peak, f.t90))
    {
        f.reason = "передний фронт не прослеживается: кривая начинается выше "
                   "10 % пика или нарастание не пересекает уровни 10/50/90 %";
        return f;
    }
    f.rise = f.t90 - f.t10;
    if (!(f.rise > 0.0))
    {
        f.reason = "длительность фронта 10-90 % не положительна "
                   "(уровни пересечены не по порядку)";
        return f;
    }
    // Спад может не попасть в окно записи — это не повод отказываться
    // от сравнения амплитуды и фронта, поэтому признак необязателен.
    double t_fall = 0.0;
    if (fall_crossing(c, i_peak, kLevelFall * f.peak, t_fall))
    {
        f.fall       = t_fall - f.t_peak;
        f.fall_valid = (f.fall > 0.0);
    }
    f.valid = true;
    return f;
}

// ГРУБЫЕ ворота. Тест обязан падать только на том, что означает настоящую
// поломку (кривая вдвое не та), а не на количественном несовпадении с
// приборным профилем.
//
// rise_is_gate — отдельный флаг, потому что применимость ворот по фронту
// зависит от КАЧЕСТВА ОЦИФРОВКИ эталона, а не от физики. У BC6 весь передний
// фронт эталона уместился в один сегмент оцифровки шириной ~9 нс, что на
// порядок-два меньше размазывания фронта ядром SPH; отношение вышло бы
// заведомо за любые разумные ворота, не сообщая ничего о расчёте. Функция
// rise_gate_ok() при этом ВСЕГДА возвращает честный числовой результат —
// флаг читает вызывающая сторона, и самопроверка остаётся содержательной.
struct ShapeGates
{
    double amp_lo  = 0.5;
    double amp_hi  = 2.0;
    double rise_lo = 1.0 / 3.0;
    double rise_hi = 3.0;
    bool   rise_is_gate = true;
};

struct ShapeMatch
{
    bool          comparable = false;
    std::string   reason;
    CurveFeatures calc;
    CurveFeatures ref;
    // Справочная величина: на сколько сдвинут расчёт при совмещении по t50.
    // НЕ ошибка расчёта — начало отсчёта эталона произвольно.
    double shift            = 0.0;   // t50_ref - t50_calc, с
    double amp_ratio        = 0.0;   // расчёт / эталон
    double rise_ratio       = 0.0;   // расчёт / эталон
    double fall_ratio       = 0.0;   // расчёт / эталон
    bool   fall_ratio_valid = false;
    double peak_offset      = 0.0;   // (t_peak_calc + shift) - t_peak_ref, с
};

inline ShapeMatch compare_curves(const Curve& calc, const Curve& ref)
{
    ShapeMatch m;
    m.calc = extract_features(calc);
    m.ref  = extract_features(ref);
    if (!m.ref.valid)
    {
        m.reason = "эталон: " + m.ref.reason;
        return m;
    }
    if (!m.calc.valid)
    {
        m.reason = "расчёт: " + m.calc.reason;
        return m;
    }
    m.shift       = m.ref.t50 - m.calc.t50;
    m.amp_ratio   = m.calc.peak / m.ref.peak;
    m.rise_ratio  = m.calc.rise / m.ref.rise;
    m.peak_offset = (m.calc.t_peak + m.shift) - m.ref.t_peak;
    if (m.calc.fall_valid && m.ref.fall_valid)
    {
        m.fall_ratio       = m.calc.fall / m.ref.fall;
        m.fall_ratio_valid = true;
    }
    m.comparable = true;
    return m;
}

inline bool amp_gate_ok(const ShapeMatch& m, const ShapeGates& g)
{
    return m.comparable && m.amp_ratio >= g.amp_lo && m.amp_ratio <= g.amp_hi;
}

inline bool rise_gate_ok(const ShapeMatch& m, const ShapeGates& g)
{
    return m.comparable && m.rise_ratio >= g.rise_lo && m.rise_ratio <= g.rise_hi;
}

// Метрики печатаются ВСЕГДА, независимо от вердикта: без них расхождение
// не диагностируется, а «зелёный» прогон не проверить глазами.
inline void print_features(const char* who, const CurveFeatures& f)
{
    constexpr double us = 1e6;
    if (!f.valid)
    {
        std::printf("  info : признаки %s: не извлечены (%s)\n", who, f.reason.c_str());
        return;
    }
    std::printf("  info : признаки %s: пик %.1f м/с при %.4f мкс; "
                "t10/t50/t90 = %.4f/%.4f/%.4f мкс; фронт 10-90 %% = %.4f мкс; ",
                who, f.peak, f.t_peak * us, f.t10 * us, f.t50 * us, f.t90 * us,
                f.rise * us);
    if (f.fall_valid) std::printf("спад до 70 %% = %.4f мкс\n", f.fall * us);
    else              std::printf("спад до 70 %% = н/д (не попал в окно записи)\n");
}

// shift_note — текст оговорки про справочный характер сдвига. Он кейсовый:
// цитата источника эталона у BC2 и BC6 разная.
inline void print_shape_report(const ShapeMatch& m, const char* shift_note)
{
    constexpr double us = 1e6;
    print_features("эталона", m.ref);
    print_features("расчёта", m.calc);
    if (!m.comparable)
    {
        std::printf("  info : сравнение невозможно: %s\n", m.reason.c_str());
        return;
    }

    std::printf("  info : совмещение по t50 переднего фронта: расчёт сдвинут на %+.4f мкс\n",
                m.shift * us);
    if (shift_note && *shift_note) std::printf("%s", shift_note);
    std::printf("  info : отношения расчёт/эталон: амплитуда = %.3f, фронт 10-90 %% = %.3f, ",
                m.amp_ratio, m.rise_ratio);
    if (m.fall_ratio_valid) std::printf("спад до 70 %% = %.3f\n", m.fall_ratio);
    else                    std::printf("спад до 70 %% = н/д\n");
    std::printf("  info : смещение пика после совмещения = %+.4f мкс\n", m.peak_offset * us);
}

// ===========================================================================
// Самопроверка компаратора
// ===========================================================================
// Требование отдельное от «тест не падает»: компаратор — измерительный
// инструмент, и он сам обязан быть проверен. Проверки работают ТОЛЬКО
// с данными эталона (прогон расчёта не нужен), поэтому выполняются в обоих
// режимах и стоят почти ничего.
//
// Проверяется ровно то, на чём стоит схема сверки:
//   * тождество            — кривая против самой себя даёт отношения 1,0 и
//                            нулевое смещение;
//   * чувствительность к амплитуде — множитель 3 обязан провалить ворота;
//   * НЕчувствительность к сдвигу времени — сдвиг на 5 мкс не меняет ничего
//     (суть схемы: начало отсчёта эталона произвольно);
//   * чувствительность к растяжению — растяжение времени в 5 раз обязано
//     провалить ворота по фронту.
//
// Последняя проверка идёт по ЧИСЛОВЫМ воротам rise_gate_ok независимо от
// gates.rise_is_gate: выключение ворот — решение о вердикте теста, а не о
// работоспособности метрики, и проверять надо именно метрику.
inline void comparator_selfcheck(const Curve& ref, const ShapeGates& gates)
{
    std::printf("\n[компаратор] Самопроверка на данных эталона\n");

    if (ref.t.size() < 2)
    {
        check(false, "самопроверка компаратора: эталон пуст или слишком короткий");
        return;
    }

    const CurveFeatures f = extract_features(ref);
    if (!check(f.valid, "признаки эталона извлекаются"))
    {
        std::printf("         причина: %s\n", f.reason.c_str());
        return;
    }
    print_features("эталона", f);

    // Допуски самопроверки — чисто арифметические: сдвиг и растяжение времени
    // выполняются в double, поэтому точное равенство ожидать нельзя нигде,
    // кроме тождества (там входные данные побитово те же).
    constexpr double kTolRatio = 1e-9;    // безразмерные отношения
    constexpr double kTolTime  = 1e-12;   // с (1e-6 мкс)

    // 1. Тождество.
    {
        const ShapeMatch m = compare_curves(ref, ref);
        check(m.comparable && m.amp_ratio == 1.0 && m.rise_ratio == 1.0 &&
                  m.fall_ratio_valid && m.fall_ratio == 1.0 &&
                  m.shift == 0.0 && m.peak_offset == 0.0,
              "тождество: эталон против самого себя даёт отношения ровно 1,0 "
              "и нулевое смещение");
        if (m.comparable)
        {
            std::printf("         амплитуда = %.17g, фронт = %.17g, спад = %.17g, "
                        "сдвиг = %.17g с, смещение пика = %.17g с\n",
                        m.amp_ratio, m.rise_ratio, m.fall_ratio, m.shift, m.peak_offset);
        }
    }

    // 2. Чувствительность к амплитуде: x3 обязано провалить ворота.
    {
        Curve scaled = ref;
        for (double& v : scaled.v) v *= 3.0;
        const ShapeMatch m = compare_curves(scaled, ref);
        check(m.comparable && !amp_gate_ok(m, gates),
              "чувствительность к амплитуде: кривая x3 проваливает ворота по амплитуде");
        // Масштаб по величине не должен затрагивать временные признаки:
        // уровни берутся в долях пика, поэтому фронт обязан остаться прежним.
        check(m.comparable && rise_gate_ok(m, gates) &&
                  std::fabs(m.rise_ratio - 1.0) <= kTolRatio,
              "избирательность: множитель по амплитуде не меняет длительность фронта");
        std::printf("         амплитуда = %.6f (ворота [%.3f, %.3f]), фронт = %.6f\n",
                    m.amp_ratio, gates.amp_lo, gates.amp_hi, m.rise_ratio);
    }

    // 3. НЕчувствительность к сдвигу времени — суть схемы сверки.
    {
        constexpr double kShift = 5e-6;   // произвольные 5 мкс
        Curve shifted = ref;
        for (double& t : shifted.t) t += kShift;
        const ShapeMatch m = compare_curves(shifted, ref);
        check(m.comparable &&
                  std::fabs(m.amp_ratio - 1.0) <= kTolRatio &&
                  std::fabs(m.rise_ratio - 1.0) <= kTolRatio &&
                  m.fall_ratio_valid && std::fabs(m.fall_ratio - 1.0) <= kTolRatio &&
                  std::fabs(m.peak_offset) <= kTolTime &&
                  amp_gate_ok(m, gates) && rise_gate_ok(m, gates),
              "нечувствительность к сдвигу времени: кривая, сдвинутая на 5 мкс, "
              "даёт те же отношения");
        if (m.comparable)
        {
            std::printf("         амплитуда = %.12f, фронт = %.12f, спад = %.12f, "
                        "смещение пика = %.3g с; справочный сдвиг = %.6f мкс\n",
                        m.amp_ratio, m.rise_ratio, m.fall_ratio, m.peak_offset,
                        m.shift * 1e6);
        }
    }

    // 4. Чувствительность к растяжению времени: x5 обязано провалить фронт.
    {
        constexpr double kStretch = 5.0;
        Curve stretched = ref;
        const double t0 = ref.t.front();
        for (double& t : stretched.t) t = t0 + kStretch * (t - t0);
        const ShapeMatch m = compare_curves(stretched, ref);
        check(m.comparable && !rise_gate_ok(m, gates),
              "чувствительность к растяжению: кривая, растянутая по времени в 5 раз, "
              "проваливает числовые ворота по фронту");
        // Растяжение по времени не трогает величины, значит амплитудные ворота
        // обязаны остаться пройденными: ворота различают род дефекта.
        check(m.comparable && amp_gate_ok(m, gates) &&
                  std::fabs(m.amp_ratio - 1.0) <= kTolRatio,
              "избирательность: растяжение по времени не меняет отношение амплитуд");
        std::printf("         фронт = %.6f (ворота [%.3f, %.3f]), амплитуда = %.6f\n",
                    m.rise_ratio, gates.rise_lo, gates.rise_hi, m.amp_ratio);
    }
}

} // namespace e2e
