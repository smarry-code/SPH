// bc6_shot.cpp
// Экспериментальная утилита кейса BC6 -- откольного опыта на карбиде бора.
// Постановка -- cases/bc6/experiment/ (копия regression/, см. cases/README.md), драйвер --
// cases/common/ceramic_shot.hpp (ключи, каталоги, чекпоинты, shot.txt, runs.csv).
//
// Сверка с приборным эталоном (cases/bc6/experiment/reference.txt), как в tests/test_e2e_bc6.cpp:
//   * лесенка расчёта (e2e_staircase.hpp): приход фронта по t50 против bc6::t_front_expected,
//     плато 1, отношения ступеней, периоды, возвраты -- со значениями допусков из bc6_case.hpp;
//   * лесенка эталона теми же средствами -- для сравнения уровней;
//   * форма (e2e_shape.hpp): отношение амплитуд, совмещение по фронту.
// Ворот нет: значения и признак «в допуске» пишутся в shot.txt.
//
// Пример:
//   OMP_NUM_THREADS=4 ./build/bc6_shot --contact-type=3 --threads=4 --out=runs/bc6

#include "ceramic_shot.hpp"
#include "e2e_shape.hpp"
#include "e2e_staircase.hpp"
#include "bc6_case.hpp"
#include "bc6_geometry.hpp"

static_assert(dim == 2, "bc6_shot требует сборки с -DSPH_DIM=2");
static_assert(axis_sim, "bc6_shot требует axis_sim = true (осесимметричная постановка)");

#ifndef SPH_CASES_DIR
#define SPH_CASES_DIR "cases"
#endif

namespace
{

using ceramic::fmtg;
using ceramic::KV;

const char* kBodyName[4] = {"foam", "impactor", "sample", "window"};

ceramic::Built build(ParticleSystem& sys, double dx)
{
    const bc6::Geometry g = bc6::build_geometry(sys, dx);
    ceramic::Built b;
    b.ntotal = g.ntotal;
    b.moni = g.moni_axis;
    for (size_t k = 0; k < 4; ++k)
    {
        const std::string p = std::string(kBodyName[k]) + "_";
        b.kv.emplace_back(p + "count", std::to_string(g.count[k]));
        b.kv.emplace_back(p + "nr_nz", std::to_string(g.n_r[k]) + "x" + std::to_string(g.n_z[k]));
    }
    b.kv.emplace_back("impact_gap", fmtg(g.impact_gap, 6));
    b.kv.emplace_back("moni_ctl_particle", std::to_string(g.moni_ctl));
    b.kv.emplace_back("moni_ctl_r", fmtg(g.moni_ctl_r, 6));
    b.kv.emplace_back("t_front_expected", fmtg(bc6::t_front_expected(g), 6));
    return b;
}

void staircase_kv(const char* who, const e2e::StaircaseFeatures& f, KV& kv)
{
    const std::string p = std::string(who) + "_";
    kv.emplace_back(p + "stair_valid", f.valid ? "1" : "0");
    if (!f.valid)
    {
        kv.emplace_back(p + "stair_reason", f.reason);
        return;
    }
    kv.emplace_back(p + "peak", fmtg(f.peak, 6));
    kv.emplace_back(p + "t_front10", fmtg(f.t_front, 6));
    kv.emplace_back(p + "t_front50", fmtg(f.t_mid, 6));
    kv.emplace_back(p + "dips", std::to_string(f.dips));
    kv.emplace_back(p + "plateaus", std::to_string(f.plateaus));
    for (size_t k = 0; k < f.plateaus && k < 3; ++k)
        kv.emplace_back(p + "plateau" + std::to_string(k + 1), fmtg(f.plateau[k], 6));
    for (size_t k = 0; k < f.dips && k < 2; ++k)
    {
        const std::string q = p + "dip" + std::to_string(k + 1) + "_";
        kv.emplace_back(q + "t", fmtg(f.t_dip[k], 6));
        kv.emplace_back(q + "v", fmtg(f.v_dip[k], 6));
        kv.emplace_back(p + "period" + std::to_string(k + 1), fmtg(f.period[k], 6));
    }
    for (size_t k = 0; k + 1 < f.plateaus && k < 2; ++k)
        kv.emplace_back(p + "ratio" + std::to_string(k + 1) + std::to_string(k + 2), fmtg(f.ratio[k], 6));
}

void analyze(const e2e::Curve& calc, const e2e::Curve& ref, double dx, size_t ntotal, KV& kv, e2e::PlotLabels& lab)
{
    lab.title = "Осевая скорость частицы-монитора, кейс BC6";
    lab.aria = "Осевая скорость частицы-монитора BC6: расчёт и приборный эталон";
    lab.grid_step = 50.0;
    lab.x_tick = 0.5;

    // Ожидаемый приход фронта -- по геометрии без частиц.
    const bc6::Geometry g = bc6::plan_geometry(dx);
    const double t_exp = bc6::t_front_expected(g);

    const e2e::StaircaseFeatures fc = e2e::extract_staircase(calc);
    e2e::print_staircase("расчёта", fc);
    staircase_kv("calc", fc, kv);
    if (!ref.empty())
    {
        const e2e::StaircaseFeatures fr = e2e::extract_staircase(ref);
        staircase_kv("ref", fr, kv);
    }

    // Допуски постановки (§6 notes/contact/06-e2e-bc6.md) -- справочно, как «в допуске: 1/0».
    if (fc.valid)
    {
        kv.emplace_back("front_ok", std::fabs(fc.t_mid - t_exp) <= bc6::kFrontTol ? "1" : "0");
        if (fc.plateaus >= 1)
            kv.emplace_back("plateau1_ok",
                            fc.plateau[0] >= bc6::kPlateau1Lo && fc.plateau[0] <= bc6::kPlateau1Hi ? "1" : "0");
        for (size_t k = 0; k + 1 < fc.plateaus && k < 2; ++k)
            kv.emplace_back("ratio" + std::to_string(k + 1) + std::to_string(k + 2) + "_ok",
                            std::fabs(fc.ratio[k] - bc6::kStepRatio) <= bc6::kStepRatioTol ? "1" : "0");
        for (size_t k = 0; k < fc.dips && k < 2; ++k)
            kv.emplace_back("period" + std::to_string(k + 1) + "_ok",
                            std::fabs(fc.period[k] - bc6::kPeriod) <= bc6::kPeriodTol ? "1" : "0");
    }

    const e2e::ShapeMatch m = e2e::compare_curves(calc, ref);
    kv.emplace_back("shape_comparable", m.comparable ? "1" : "0");
    if (m.comparable)
    {
        kv.emplace_back("amp_ratio", fmtg(m.amp_ratio, 6));
        kv.emplace_back("rise_ratio", fmtg(m.rise_ratio, 6));
        kv.emplace_back("fall_ratio", m.fall_ratio_valid ? fmtg(m.fall_ratio, 6) : "nan");
        kv.emplace_back("amp_in_test_gate", m.amp_ratio >= 0.5 && m.amp_ratio <= 2.0 ? "1" : "0");
        // Сдвиг эталона для графика -- по t10, как reference_shift() в tests/test_e2e_bc6.cpp.
        if (m.calc.valid && m.ref.valid) lab.ref_shift = m.calc.t10 - m.ref.t10;
        kv.emplace_back("ref_shift_t10", fmtg(lab.ref_shift, 6));
    }
    else
    {
        kv.emplace_back("shape_reason", m.reason);
    }

    std::ostringstream os;
    os << "Граница «керамика / окно», ось. Расчёт: delta_x = " << e2e::ru_num(dx * 1e3, 2) << " мм, "
       << e2e::ru_num(static_cast<double>(ntotal), 0) << " "
       << e2e::ru_plural(ntotal, "частица", "частицы", "частиц") << ". ";
    if (lab.ref_shift != 0.0) os << "Эталон совмещён по фронту (" << e2e::ru_num(lab.ref_shift * 1e6, 2) << " мкс). ";
    os << "Крестик -- прогон оборван до t_total.";
    lab.subtitle = os.str();
}

} // namespace

int main(int argc, char** argv)
{
    ceramic::CaseHooks h;
    h.name = "bc6";
    h.title = "BC6: откольный опыт на карбиде бора";
    h.setup_dir = std::filesystem::path(SPH_CASES_DIR) / "bc6" / "experiment";
    h.dx_default = bc6::kDxFull;
    h.build = build;
    h.analyze = analyze;
    return ceramic::run_main(argc, argv, h);
}
