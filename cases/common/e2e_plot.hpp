#pragma once
// e2e_plot.hpp
// Сводный SVG-график кривых монитора: эталон и расчёт по типам контакта.
//
// Происхождение — tests/test_e2e_bc2.cpp, раздел 8 (emit_plot).
// Сам BC2 на этот заголовок не переключён, см. пояснение в e2e_harness.hpp.
//
// Что параметризовано против оригинала:
//   * эталон передаётся аргументом (BC2 читал его сам внутри emit_plot);
//   * геометрия кейса больше не нужна — подпись собирает вызывающая сторона;
//   * заголовок, подпись, aria-label, названия осей, шаг сетки по Y и по X,
//     запасы по осям — поля PlotLabels (были литералами);
//   * ref_shift — сдвиг эталона по времени при отрисовке (у BC2 не было).
//
// ПРО ref_shift. Начало отсчёта времени у приборного эталона произвольно, и
// это относится к обоим кейсам. BC2 рисует эталон на сырой абсолютной оси и
// объясняет видимый разрыв текстом отчёта. Для BC6 постановка (§5.3) прямо
// требует вести сверку только после выравнивания по фронту, поэтому у него
// эталон уместнее рисовать совмещённым: иначе картинка показывает разрыв,
// который сама же постановка объявила лишённым смысла. Величина сдвига
// обязана быть названа в подписи — иначе график вводит в заблуждение уже в
// другую сторону.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "e2e_curve.hpp"
#include "e2e_harness.hpp"

namespace e2e
{

struct PlotLabels
{
    std::string title;        // заголовок над полем
    std::string subtitle;     // строка пояснения под заголовком
    std::string aria;         // aria-label для доступности
    std::string ref_legend = "эталон (прибор)";
    std::string x_axis     = "время, мкс";
    std::string y_axis     = "осевая скорость, м/с";
    double grid_step  = 300.0;   // шаг линий сетки по Y, м/с
    double x_tick     = 0.5;     // шаг подписей по X, мкс
    double y_headroom = 1.08;
    double x_headroom = 1.02;
    double ref_shift  = 0.0;     // с, прибавляется ко времени эталона при отрисовке
};

inline void emit_plot(const std::vector<PlotSeries>& series, const Curve& reference,
                      const PlotLabels& lab, const std::filesystem::path& out_path)
{
    namespace fs = std::filesystem;
    std::printf("\n[график] %s\n", out_path.string().c_str());

    bool any_curve = false;
    for (const PlotSeries& s : series) any_curve = any_curve || !s.curve.empty();
    if (!any_curve)
    {
        std::printf("  info : кривых расчёта нет, график не построен\n");
        return;
    }

    Curve ref = reference;
    if (lab.ref_shift != 0.0)
    {
        for (double& t : ref.t) t += lab.ref_shift;
    }

    std::vector<PlotSeries> plot = series;
    for (PlotSeries& s : plot) s.curve = thin_curve(s.curve);

    // Стиль серии. Одного цвета мало: типы 0 и 1 на некоторых кейсах совпадают
    // до четвёртой значащей цифры, и сплошная линия, нарисованная позже,
    // закрывает раннюю целиком — легенда обещала бы четыре кривые, а на поле
    // их было бы три. Поэтому к цвету добавлены штрих и толщина: нижняя
    // сплошная просматривается сквозь пунктир верхней. Побочно это даёт
    // второй признак помимо цвета, что нужно и при цветовой слепоте.
    // Ключ — ТИП, а не место в списке: если какого-то типа нет (прогон не
    // делали), остальные не должны сменить штрих.
    struct SeriesStyle { double width; const char* dash; };
    const auto style_of = [](int type) -> SeriesStyle
    {
        switch (type & 3)
        {
            case 0:  return {3.0, ""};
            case 1:  return {2.4, "9 6"};
            case 2:  return {2.0, ""};
            default: return {2.2, "2.5 4"};
        }
    };

    // Пределы осей. Потолок по Y — ближайшее сверху кратное шагу сетки, чтобы
    // линии попадали на круглые значения и кривая не упиралась в рамку.
    double vmax = 0.0, tmax = 0.0;
    for (const PlotSeries& s : plot)
    {
        for (size_t i = 0; i < s.curve.t.size(); ++i)
        {
            vmax = std::max(vmax, s.curve.v[i]);
            tmax = std::max(tmax, s.curve.t[i]);
        }
    }
    for (size_t i = 0; i < ref.t.size(); ++i)
    {
        vmax = std::max(vmax, ref.v[i]);
        tmax = std::max(tmax, ref.t[i]);
    }
    const double kGrid = (lab.grid_step > 0.0) ? lab.grid_step : 300.0;
    const double ymax = std::ceil(std::max(vmax, 1.0) * lab.y_headroom / kGrid) * kGrid;
    const double xmax_us = std::max(tmax * 1e6, 1e-9) * lab.x_headroom;

    // Поле справа расширено под прямые подписи кривых: в светлой теме аква и
    // жёлтый идут ниже 3:1 к фону, и подпись у конца кривой — обязательная
    // страховка, а не украшение. Сверху добавлено место под строку легенды.
    constexpr double W = 880.0, H = 470.0;
    constexpr double ML = 68.0, MR = 104.0, MT = 84.0, MB = 56.0;
    const double PW = W - ML - MR, PH = H - MT - MB;

    const auto sx = [&](double t_us) { return ML + t_us / xmax_us * PW; };
    const auto sy = [&](double v) { return MT + PH - std::min(std::max(v, 0.0), ymax) / ymax * PH; };

    std::ostringstream os;
    os.setf(std::ios::fixed);
    const auto path_of = [&](const Curve& c)
    {
        std::ostringstream d;
        d.setf(std::ios::fixed);
        d.precision(1);
        for (size_t i = 0; i < c.t.size(); ++i)
        {
            d << (i == 0 ? 'M' : 'L') << sx(c.t[i] * 1e6) << ',' << sy(c.v[i]) << ' ';
        }
        return d.str();
    };

    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 880 470\" width=\"880\" height=\"470\" "
          "role=\"img\" aria-label=\"" << lab.aria << "\">\n";
    // Палитра задаётся токенами на :root-уровне svg и переопределяется для
    // тёмной темы, поэтому график читается в обеих.
    os << "<style>\n"
          "svg { --surface-1:#fcfcfb; --text-primary:#0b0b0b; --text-secondary:#52514e;\n"
          "      --grid:#e4e3df; --axis:#b9b8b2; --ref:#0b0b0b;\n"
          "      --s0:#2a78d6; --s1:#eb6834; --s2:#1baf7a; --s3:#eda100;\n"
          "      font-family:-apple-system,BlinkMacSystemFont,\"Segoe UI\",Roboto,sans-serif; }\n"
          "@media (prefers-color-scheme: dark) {\n"
          "  svg { --surface-1:#1a1a19; --text-primary:#ffffff; --text-secondary:#c3c2b7;\n"
          "        --grid:#333330; --axis:#5a5955; --ref:#ffffff;\n"
          "        --s0:#3987e5; --s1:#d95926; --s2:#199e70; --s3:#c98500; }\n"
          "}\n"
          ".ttl { fill:var(--text-primary); font-size:15px; font-weight:600; }\n"
          ".sub { fill:var(--text-secondary); font-size:11.5px; }\n"
          ".axl { fill:var(--text-secondary); font-size:11px; }\n"
          ".tik { fill:var(--text-secondary); font-size:10.5px; }\n"
          ".lbl { fill:var(--text-primary); font-size:11px; font-weight:600; }\n"
          ".leg { fill:var(--text-secondary); font-size:11px; }\n"
          ".ln  { fill:none; stroke-width:2; stroke-linejoin:round; stroke-linecap:round; }\n"
          "</style>\n";
    os.precision(1);
    os << "<rect width=\"880\" height=\"470\" fill=\"var(--surface-1)\"/>\n";
    os << "<text class=\"ttl\" x=\"" << ML << "\" y=\"24\">" << lab.title << "</text>\n";
    os << "<text class=\"sub\" x=\"" << ML << "\" y=\"41\">" << lab.subtitle << "</text>\n";

    // Легенда — строкой над полем, колонки с постоянным шагом. Текст носит
    // текстовые токены, цвет несёт только штрих слева от подписи, поэтому
    // идентичность не держится на одном цвете.
    {
        constexpr double kLegY = 62.0, kLegPitch = 132.0, kLegRule = 22.0, kLegGap = 6.0;
        // Штрих и толщина в легенде — ровно те же, что у кривой на поле:
        // иначе легенда опознаёт серию не по тому признаку, по которому её
        // приходится узнавать на графике.
        const auto legend_entry = [&](double lx, const std::string& stroke,
                                      const SeriesStyle& st, bool faded,
                                      const std::string& text)
        {
            os << "<line x1=\"" << lx << "\" y1=\"" << kLegY
               << "\" x2=\"" << (lx + kLegRule) << "\" y2=\"" << kLegY
               << "\" stroke=\"" << stroke << "\" stroke-width=\"" << st.width << "\"";
            if (st.dash[0] != '\0') os << " stroke-dasharray=\"" << st.dash << "\"";
            if (faded) os << " opacity=\"0.85\"";
            os << "/>\n";
            os << "<text class=\"leg\" x=\"" << (lx + kLegRule + kLegGap)
               << "\" y=\"" << (kLegY + 3.8) << "\">" << text << "</text>\n";
        };
        double lx = ML;
        for (size_t i = 0; i < plot.size(); ++i)
        {
            std::ostringstream t;
            t << "тип " << plot[i].type << " — " << contact_type_name(plot[i].type);
            legend_entry(lx, "var(--s" + std::to_string(plot[i].type) + ")",
                         style_of(plot[i].type), false, t.str());
            lx += kLegPitch;
        }
        if (!ref.empty()) legend_entry(lx, "var(--ref)", SeriesStyle{2.0, "6 4"}, true,
                                       lab.ref_legend);
    }

    for (double v = 0.0; v <= ymax + 1e-9; v += kGrid)
    {
        const double y = sy(v);
        os << "<line x1=\"" << ML << "\" y1=\"" << y << "\" x2=\"" << (ML + PW) << "\" y2=\"" << y
           << "\" stroke=\"var(--grid)\" stroke-width=\"1\"/>\n";
        os << "<text class=\"tik\" x=\"" << (ML - 10) << "\" y=\"" << (y + 3.6)
           << "\" text-anchor=\"end\">" << ru_num(v, 0) << "</text>\n";
    }
    const double x_tick = (lab.x_tick > 0.0) ? lab.x_tick : 0.5;
    for (double t = 0.0; t <= xmax_us + 1e-9; t += x_tick)
    {
        os << "<text class=\"tik\" x=\"" << sx(t) << "\" y=\"" << (MT + PH + 18)
           << "\" text-anchor=\"middle\">" << ru_num(t, 1) << "</text>\n";
    }
    os << "<line x1=\"" << ML << "\" y1=\"" << (MT + PH) << "\" x2=\"" << (ML + PW) << "\" y2=\"" << (MT + PH)
       << "\" stroke=\"var(--axis)\" stroke-width=\"1\"/>\n";
    os << "<text class=\"axl\" x=\"" << (ML + PW / 2) << "\" y=\"" << (H - 14)
       << "\" text-anchor=\"middle\">" << lab.x_axis << "</text>\n";
    os << "<text class=\"axl\" x=\"18\" y=\"" << (MT + PH / 2) << "\" text-anchor=\"middle\" "
       << "transform=\"rotate(-90 18 " << (MT + PH / 2) << ")\">" << lab.y_axis << "</text>\n";

    if (!ref.empty())
    {
        os << "<path class=\"ln\" d=\"" << path_of(ref)
           << "\" stroke=\"var(--ref)\" stroke-dasharray=\"6 4\" opacity=\"0.85\"/>\n";
    }
    // Кривые расчёта: стиль — из style_of, тот же, что в легенде.
    for (size_t i = 0; i < plot.size(); ++i)
    {
        if (plot[i].curve.empty()) continue;
        const SeriesStyle st = style_of(plot[i].type);
        os << "<path class=\"ln\" d=\"" << path_of(plot[i].curve)
           << "\" stroke=\"var(--s" << plot[i].type << ")\" stroke-width=\"" << st.width << "\"";
        if (st.dash[0] != '\0') os << " stroke-dasharray=\"" << st.dash << "\"";
        os << "/>\n";
    }

    // Прямые подписи у концов кривых с разведением по вертикали: цвет в
    // светлой теме ниже 3:1 у двух слотов, и подпись — обязательная страховка.
    struct EndLabel
    {
        double y_end;     // конец кривой, куда ставится крестик обрыва
        double y_label;   // та же точка после разведения подписей
        double x_end;
        int    type;
        bool   aborted;
    };
    std::vector<EndLabel> ends;
    for (const PlotSeries& s : plot)
    {
        if (s.curve.empty()) continue;
        const double ye = sy(s.curve.v.back());
        ends.push_back({ye, ye, sx(s.curve.t.back() * 1e6), s.type, s.aborted});
    }
    std::sort(ends.begin(), ends.end(),
              [](const EndLabel& a, const EndLabel& b) { return a.y_label < b.y_label; });
    // Разведение по вертикали в два прохода. Один проход вниз недостаточен:
    // когда кривые кончаются рядом у нижнего края (smoke-режим — все четыре
    // у нуля), подписи уезжают за рамку на подписи оси времени. Обратный
    // проход возвращает пачку внутрь поля.
    constexpr double kLabelMin = 14.0;
    const double y_lo = MT + 6.0, y_hi = MT + PH - 6.0;
    for (size_t i = 1; i < ends.size(); ++i)
    {
        if (ends[i].y_label - ends[i - 1].y_label < kLabelMin)
        {
            ends[i].y_label = ends[i - 1].y_label + kLabelMin;
        }
    }
    if (!ends.empty() && ends.back().y_label > y_hi)
    {
        ends.back().y_label = y_hi;
        for (size_t i = ends.size() - 1; i > 0; --i)
        {
            if (ends[i].y_label - ends[i - 1].y_label < kLabelMin)
            {
                ends[i - 1].y_label = ends[i].y_label - kLabelMin;
            }
        }
        for (EndLabel& e : ends) e.y_label = std::max(e.y_label, y_lo);
    }
    for (const EndLabel& e : ends)
    {
        // Крестик в точке обрыва: обрыв — свойство данных и обязан быть виден
        // на самом графике, а не только в подписи сбоку.
        if (e.aborted)
        {
            os << "<path d=\"M" << (e.x_end - 4) << "," << (e.y_end - 4)
               << "L" << (e.x_end + 4) << "," << (e.y_end + 4)
               << "M" << (e.x_end + 4) << "," << (e.y_end - 4)
               << "L" << (e.x_end - 4) << "," << (e.y_end + 4)
               << "\" stroke=\"var(--s" << e.type << ")\" stroke-width=\"2\" fill=\"none\"/>\n";
        }
        os << "<line x1=\"" << (ML + PW + 6) << "\" y1=\"" << e.y_label
           << "\" x2=\"" << (ML + PW + 24) << "\" y2=\"" << e.y_label
           << "\" stroke=\"var(--s" << e.type << ")\" stroke-width=\"2.5\"/>\n";
        os << "<text class=\"lbl\" x=\"" << (ML + PW + 30) << "\" y=\"" << (e.y_label + 3.8)
           << "\">тип " << e.type << (e.aborted ? " ✕" : "") << "</text>\n";
    }
    os << "</svg>\n";

    std::error_code ec;
    if (out_path.has_parent_path()) fs::create_directories(out_path.parent_path(), ec);
    std::ofstream f(out_path, std::ios::binary | std::ios::trunc);
    if (!f)
    {
        std::printf("  FAIL : не удалось открыть файл для записи графика\n");
        ++failures;
        return;
    }
    const std::string svg = os.str();
    f << svg;
    f.close();
    if (!f)
    {
        std::printf("  FAIL : ошибка записи графика\n");
        ++failures;
        return;
    }
    size_t pts = 0;
    for (const PlotSeries& s : plot) pts += s.curve.t.size();
    std::printf("  ok   : график записан (%zu %s, %zu точек расчёта, %zu эталона, %zu байт)\n",
                plot.size(), ru_plural(plot.size(), "кривая", "кривые", "кривых"),
                pts, ref.t.size(), svg.size());
}

} // namespace e2e
