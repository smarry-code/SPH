// bc2_shot.cpp
// Экспериментальная утилита кейса BC2 -- удар по карбиду бора (Grady & Moody, SAND96-0551).
// Постановка -- cases/bc2/experiment/ (копия regression/, см. cases/README.md), драйвер --
// cases/common/ceramic_shot.hpp (ключи, каталоги, чекпоинты, shot.txt, runs.csv).
//
// Сверка с приборным эталоном (cases/bc2/experiment/reference.txt) -- по форме кривой, как в
// tests/test_e2e_bc2.cpp: признаки фронта и пика, совмещение по t50 (начало отсчёта эталона
// произвольно), отношения амплитуды, фронта 10-90 % и спада. Ворот нет: значения пишутся в shot.txt.
//
// Пример:
//   OMP_NUM_THREADS=4 ./build/bc2_shot --contact-type=3 --threads=4 --out=runs/bc2

#include "ceramic_shot.hpp"
#include "e2e_shape.hpp"
#include "bc2_geometry.hpp"

static_assert(dim == 2, "bc2_shot требует сборки с -DSPH_DIM=2");
static_assert(axis_sim, "bc2_shot требует axis_sim = true (осесимметричная постановка)");

#ifndef SPH_CASES_DIR
#define SPH_CASES_DIR "cases"
#endif

namespace
{

using ceramic::fmtg;
using ceramic::KV;

ceramic::Built build(ParticleSystem& sys, double dx)
{
    const bc2::Bc2Geometry g = bc2::build_bc2_geometry(sys, dx);
    ceramic::Built b;
    b.ntotal = g.ntotal;
    b.moni = g.moni_particle;
    for (size_t k = 0; k < g.layer_count.size(); ++k)
    {
        const std::string p = "layer" + std::to_string(k + 1) + "_";
        b.kv.emplace_back(p + "count", std::to_string(g.layer_count[k]));
        b.kv.emplace_back(p + "nr_nz", std::to_string(g.layer_nr[k]) + "x" + std::to_string(g.layer_nz[k]));
        b.kv.emplace_back(p + "dr_dz", fmtg(g.layer_dr[k], 6) + "x" + fmtg(g.layer_dz[k], 6));
    }
    return b;
}

void analyze(const e2e::Curve& calc, const e2e::Curve& ref, double dx, size_t ntotal, KV& kv, e2e::PlotLabels& lab)
{
    lab.title = "Осевая скорость частицы-монитора, кейс BC2";
    lab.aria = "Осевая скорость частицы-монитора BC2: расчёт и приборный эталон";
    lab.grid_step = 300.0;
    lab.x_tick = 0.5;

    const e2e::ShapeMatch m = e2e::compare_curves(calc, ref);
    const auto feat = [&](const char* who, const e2e::CurveFeatures& f) {
        const std::string p = std::string(who) + "_";
        kv.emplace_back(p + "valid", f.valid ? "1" : "0");
        if (!f.valid)
        {
            kv.emplace_back(p + "reason", f.reason);
            return;
        }
        kv.emplace_back(p + "peak", fmtg(f.peak, 6));
        kv.emplace_back(p + "t_peak", fmtg(f.t_peak, 6));
        kv.emplace_back(p + "t10", fmtg(f.t10, 6));
        kv.emplace_back(p + "t50", fmtg(f.t50, 6));
        kv.emplace_back(p + "t90", fmtg(f.t90, 6));
        kv.emplace_back(p + "rise_10_90", fmtg(f.rise, 6));
        kv.emplace_back(p + "fall_to_70", f.fall_valid ? fmtg(f.fall, 6) : "nan");
    };
    feat("calc", m.calc);
    feat("ref", m.ref);
    kv.emplace_back("shape_comparable", m.comparable ? "1" : "0");
    if (!m.comparable)
    {
        kv.emplace_back("shape_reason", m.reason);
    }
    else
    {
        // shift -- справочная величина: время эталона не привязано к моменту удара.
        kv.emplace_back("shift_t50", fmtg(m.shift, 6));
        kv.emplace_back("amp_ratio", fmtg(m.amp_ratio, 6));
        kv.emplace_back("rise_ratio", fmtg(m.rise_ratio, 6));
        kv.emplace_back("fall_ratio", m.fall_ratio_valid ? fmtg(m.fall_ratio, 6) : "nan");
        kv.emplace_back("peak_offset", fmtg(m.peak_offset, 6));
        const e2e::ShapeGates g;   // ворота теста BC2 (e2e_shape.hpp) -- справочно
        kv.emplace_back("amp_in_test_gate", e2e::amp_gate_ok(m, g) ? "1" : "0");
        kv.emplace_back("rise_in_test_gate", e2e::rise_gate_ok(m, g) ? "1" : "0");
        lab.ref_shift = -m.shift;   // эталон рисуется совмещённым с расчётом по t50
    }
    e2e::print_shape_report(m, "         сдвиг справочный: начало отсчёта времени эталона произвольно\n");

    std::ostringstream os;
    os << "Граница «керамика / окно», ось. Расчёт: delta_x = " << e2e::ru_num(dx * 1e3, 2) << " мм, "
       << e2e::ru_num(static_cast<double>(ntotal), 0) << " "
       << e2e::ru_plural(ntotal, "частица", "частицы", "частиц") << ". ";
    if (lab.ref_shift != 0.0) os << "Эталон совмещён по t50 фронта (" << e2e::ru_num(lab.ref_shift * 1e6, 2) << " мкс). ";
    os << "Крестик -- прогон оборван до t_total.";
    lab.subtitle = os.str();
}

} // namespace

int main(int argc, char** argv)
{
    ceramic::CaseHooks h;
    h.name = "bc2";
    h.title = "BC2: удар по карбиду бора";
    h.setup_dir = std::filesystem::path(SPH_CASES_DIR) / "bc2" / "experiment";
    h.dx_default = 2e-4;
    h.build = build;
    h.analyze = analyze;
    return ceramic::run_main(argc, argv, h);
}
