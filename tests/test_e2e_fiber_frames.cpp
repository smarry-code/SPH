// test_e2e_fiber_frames.cpp
// Самопроверка читателя кадров нитей tests/e2e_fiber_frames.hpp.
//
// Схема: небольшая модель из двух нитей с заранее известными координатами и
// скоростями (в том числе с крошечным и огромным значением) пишется боевым
// fiber::fiber_output() в собственный временный каталог (ScopedWorkDir), затем
// читается новым читателем и сравнивается с записанным. fiber_output печатает
// числа потоком по умолчанию — 6 значащих цифр, поэтому допуск сравнения
// относительный, 1e-5 (половина единицы шестого знака даёт до 5e-6).
// Затем проверяются: порядок кадров по числовому номеру шага, поиск бусины по
// id после разрыва связи (object_id меняется, id нет), живучесть кадров, а
// также обнаружение испорченного кадра ("nan" в координате; обрезанный файл).
//
// Сборка: SPH_DIM=3 (формат x[dim] осмыслен только в 3D), источники
// ${SPH_PIPELINE_SOURCES} (fiber_output.cpp входит в конвейер).
#include "fiber/fiber_output.hpp"
#include "fiber/model.hpp"
#include "param.hpp"
#include "particle_system.hpp"

#include "e2e_fiber_frames.hpp"
#include "e2e_harness.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

static_assert(dim == 3, "test_e2e_fiber_frames требует сборки с SPH_DIM=3");

namespace
{

namespace fs = std::filesystem;

constexpr double kRelTol = 1.0e-5; // 6 значащих цифр вывода по умолчанию

// Нить из n бусин вдоль dir, шаг spacing; связи (k, k+1). Индексы глобальные.
std::size_t add_chain(fiber::Model& model, std::size_t n, const double (&origin)[3], const double (&dir)[3],
                      double spacing, int object_id)
{
    const std::size_t first = model.fiber.size();
    model.fiber.resize(first + n);
    model.topology.prev.resize(first + n, -1);
    model.topology.next.resize(first + n, -1);
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t g = first + k;
        for (std::size_t d = 0; d < 3; ++d)
        {
            model.fiber.x[d][g] = origin[d] + static_cast<double>(k) * spacing * dir[d];
            model.fiber.vx[d][g] = 0.0;
        }
        model.fiber.mass[g] = 1.0e-6;
        model.fiber.object_id[g] = object_id;
        model.fiber.material_id[g] = 0;
        if (k > 0) model.topology.prev[g] = static_cast<long long>(g - 1);
        if (k + 1 < n) model.topology.next[g] = static_cast<long long>(g + 1);
    }
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        fiber::FiberBond b;
        b.i = first + k;
        b.j = first + k + 1;
        b.l0 = spacing;
        b.l_break = 1.5 * spacing;
        model.topology.bonds.push_back(b);
    }
    return first;
}

bool close_rel(double a, double b)
{
    if (a == b) return true;
    return std::fabs(a - b) <= kRelTol * std::max(std::fabs(a), std::fabs(b));
}

bool close_vec(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return close_rel(a[0], b[0]) && close_rel(a[1], b[1]) && close_rel(a[2], b[2]);
}

std::string slurp(const fs::path& p)
{
    std::ifstream in(p);
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
}

void spit(const fs::path& p, const std::string& s)
{
    std::ofstream out(p);
    out << s;
}

} // namespace

int main()
{
    std::printf("=== e2e_fiber_frames: самопроверка читателя кадров нитей ===\n");

    const fs::path cwd = fs::current_path();
    const e2e::DirSnapshot before = e2e::snapshot_dir(cwd);
    const std::size_t fiber_files_before = e2e::list_fiber_frames(cwd).size();

    // --- сцена -------------------------------------------------------------
    fiber::Model model;
    model.solid.resize(1); // как в clamp_e2e_fiber.hpp: фиктивная частица
    model.solid.inzializ_index(0);
    model.fiber.resize(0);
    model.topology.resize(0);
    model.materials.assign(1, fiber::FiberMaterial{});

    const double oA[3] = {0.0, 0.0, 0.0}, dA[3] = {1, 0, 0};
    const double oB[3] = {0.5, -0.25, 0.125}, dB[3] = {0, 1, 0};
    const std::size_t a0 = add_chain(model, 4, oA, dA, 0.01, 0); // бусины 0..3, связи 0..2
    const std::size_t b0 = add_chain(model, 3, oB, dB, 0.02, 1); // бусины 4..6, связи 3..4
    const std::size_t nbeads = model.fiber.size();
    const std::size_t nbonds0 = model.topology.bonds.size();

    // Известные скорости; бусина 2 — крошечные/огромные значения, бусина 5 — отрицательные.
    for (std::size_t k = 0; k < nbeads; ++k)
    {
        model.fiber.vx[0][k] = 0.5 + 0.25 * static_cast<double>(k);
        model.fiber.vx[1][k] = -1.0e-3 * static_cast<double>(k);
        model.fiber.vx[2][k] = 0.0;
    }
    model.fiber.vx[0][2] = 1.2345678e-300;
    model.fiber.vx[1][2] = -9.87654321e+200;
    model.fiber.vx[2][2] = 3.0e-17;
    model.fiber.vx[0][5] = -123.456789;
    model.fiber.vx[2][5] = 7.0e-4;
    model.fiber.x[2][a0] = -2.5e-5;

    {
        e2e::ScopedWorkDir work("sph_e2e_fiber_frames", "main");
        const fs::path dir = work.dir();

        // кадр 1: исходное состояние; кадр 2: сдвиг + разрыв связи (id тот же, object_id новый);
        // кадр 3: 6-значный номер шага (проверка числового порядка и разбора имени)
        const std::size_t steps[3] = {2, 10, 123456};
        const double times[3] = {1.5e-6, 2.25e-5, 0.125};

        // Копии записанного состояния для сравнения (с учётом изменений между кадрами).
        std::vector<std::array<double, 3>> wx[3], wv[3];
        std::vector<int> woid[3];
        std::size_t wbonds[3] = {};
        auto remember = [&](int f)
        {
            for (std::size_t k = 0; k < nbeads; ++k)
            {
                wx[f].push_back({model.fiber.x[0][k], model.fiber.x[1][k], model.fiber.x[2][k]});
                wv[f].push_back({model.fiber.vx[0][k], model.fiber.vx[1][k], model.fiber.vx[2][k]});
                woid[f].push_back(model.fiber.object_id[k]);
            }
            wbonds[f] = model.topology.bonds.size();
        };

        try
        {
            remember(0);
            fiber::fiber_output(model, steps[0], times[0]);

            for (std::size_t k = 0; k < nbeads; ++k) model.fiber.x[0][k] += 1.0e-3 * (1.0 + static_cast<double>(k));
            model.topology.bonds.erase(model.topology.bonds.begin() + 1); // рвём связь (1,2)
            for (std::size_t k = 2; k <= 3; ++k) model.fiber.object_id[k] = 7; // отлетевшая часть
            remember(1);
            fiber::fiber_output(model, steps[1], times[1]);

            model.fiber.vx[0][1] = 42.0;
            remember(2);
            fiber::fiber_output(model, steps[2], times[2]);
        }
        catch (const std::exception& e)
        {
            e2e::check(false, std::string("fiber_output бросил исключение: ") + e.what());
            std::printf("\n=== e2e_fiber_frames: ПРОВАЛЕН ===\n");
            return 1;
        }

        std::printf("\n[1] Перечень кадров\n");
        const std::vector<fs::path> files = e2e::list_fiber_frames(dir);
        e2e::check(files.size() == 3, "найдено 3 кадра, найдено " + std::to_string(files.size()));
        if (files.size() == 3)
        {
            e2e::check(files[0].filename() == "fiber_output_00002.txt" &&
                           files[1].filename() == "fiber_output_00010.txt" &&
                           files[2].filename() == "fiber_output_123456.txt",
                       "порядок по числовому номеру шага: 2, 10, 123456");
        }
        // посторонние файлы в каталоге не должны попадать в перечень
        spit(dir / "fiber_output_abc.txt", "x");
        spit(dir / "sph_output_00001.txt", "x");
        spit(dir / "fiber_output_00005.dat", "x");
        e2e::check(e2e::list_fiber_frames(dir).size() == 3, "посторонние имена (abc, sph_output_*, .dat) не учитываются");
        fs::remove(dir / "fiber_output_abc.txt");
        fs::remove(dir / "sph_output_00001.txt");
        fs::remove(dir / "fiber_output_00005.dat");

        std::printf("\n[2] Чтение кадров и сравнение с записанным\n");
        for (std::size_t f = 0; f < files.size() && f < 3; ++f)
        {
            const std::string tag = "кадр " + std::to_string(steps[f]);
            e2e::FiberFrame fr;
            std::string err;
            const bool ok = e2e::read_fiber_frame(files[f], fr, &err);
            e2e::check(ok, tag + ": разбор без ошибки" + (ok ? "" : " (" + err + ")"));
            if (!ok) continue;

            e2e::check(fr.step == steps[f], tag + ": STEP из заголовка = " + std::to_string(fr.step));
            e2e::check(close_rel(fr.time, times[f]), tag + ": TIME = " + std::to_string(fr.time));
            e2e::check(fr.n_particles_header == nbeads && fr.id.size() == nbeads,
                       tag + ": N_PARTICLES = " + std::to_string(fr.n_particles_header) + ", строк " +
                           std::to_string(fr.id.size()));
            e2e::check(fr.n_bonds_header == wbonds[f] && fr.bonds.size() == wbonds[f],
                       tag + ": N_BONDS = " + std::to_string(fr.n_bonds_header) + ", строк " +
                           std::to_string(fr.bonds.size()));
            if (fr.id.size() != nbeads) continue;

            bool ids_ok = true, x_ok = true, v_ok = true, oid_ok = true, mid_ok = true;
            for (std::size_t k = 0; k < nbeads; ++k)
            {
                ids_ok = ids_ok && fr.id[k] == k;
                x_ok = x_ok && close_vec(fr.x[k], wx[f][k]);
                v_ok = v_ok && close_vec(fr.v[k], wv[f][k]);
                oid_ok = oid_ok && fr.object_id[k] == woid[f][k];
                mid_ok = mid_ok && fr.material_id[k] == 0;
            }
            e2e::check(ids_ok, tag + ": id = 0..N-1 подряд");
            e2e::check(x_ok, tag + ": координаты совпали (отн. допуск 1e-5)");
            e2e::check(v_ok, tag + ": скорости совпали (в т.ч. 1.2e-300 и -9.9e+200)");
            e2e::check(oid_ok && mid_ok, tag + ": object_id / material_id совпали");

            // Связи: длина и удлинение пересчитываются независимо по координатам из кадра.
            bool bond_ok = true;
            for (const auto& b : fr.bonds)
            {
                if (b.i >= nbeads || b.j >= nbeads) { bond_ok = false; break; }
                double l2 = 0.0;
                for (std::size_t d = 0; d < 3; ++d)
                {
                    const double dc = fr.x[b.j][d] - fr.x[b.i][d];
                    l2 += dc * dc;
                }
                const double l = std::sqrt(l2);
                // координаты кадра округлены до 6 цифр (до 5e-7 м при |x|~0.5), поэтому допуск на l мягче
                bond_ok = bond_ok && std::fabs(l - b.l) <= 1.0e-3 * b.l0 &&
                          std::fabs(b.elong - (b.l - b.l0) / b.l0) <= 1.0e-3 && close_rel(b.l_break, 1.5 * b.l0);
            }
            e2e::check(bond_ok, tag + ": длины связей и удлинения согласованы с координатами");
        }

        std::printf("\n[3] Поиск бусины по id и разрыв связи\n");
        {
            e2e::FiberFrame f0, f1;
            const bool r0 = e2e::read_fiber_frame(files[0], f0, nullptr);
            const bool r1 = e2e::read_fiber_frame(files[1], f1, nullptr);
            e2e::check(r0 && r1, "кадры 2 и 10 читаются");
            if (r0 && r1)
            {
                std::array<double, 3> x{}, v{};
                const bool found = e2e::fiber_frame_bead(f0, 2, x, v);
                e2e::check(found && close_vec(x, wx[0][2]) && close_vec(v, wv[0][2]),
                           "бусина id=2 в кадре 2: значения совпали");
                e2e::check(e2e::fiber_frame_bead(f1, 5, x, v) && close_vec(v, wv[1][5]) && close_vec(x, wx[1][5]),
                           "бусина id=5 в кадре 10: значения совпали");
                e2e::check(!e2e::fiber_frame_bead(f1, nbeads, x, v), "несуществующий id (= N) не найден");
                e2e::check(!e2e::fiber_frame_bead(f1, 100000, x, v), "несуществующий id (100000) не найден");
                e2e::check(f0.object_id[2] == 0 && f1.object_id[2] == 7,
                           "после разрыва у бусины 2 новый object_id (0 -> 7), id прежний");
                e2e::check(f1.bonds.size() + 1 == f0.bonds.size() && f0.bonds.size() == nbonds0,
                           "число связей уменьшилось на одну: " + std::to_string(f0.bonds.size()) + " -> " +
                               std::to_string(f1.bonds.size()));
                // бусина, найденная по id в обоих кадрах, сместилась на известную величину
                std::array<double, 3> xa{}, va{}, xb{}, vb{};
                const bool both = e2e::fiber_frame_bead(f0, 4, xa, va) && e2e::fiber_frame_bead(f1, 4, xb, vb);
                e2e::check(both && std::fabs((xb[0] - xa[0]) - 1.0e-3 * 5.0) <= 2.0e-6,
                           "смещение бусины id=4 между кадрами = 5e-3");
            }
        }

        std::printf("\n[4] Живучесть: исправные кадры\n");
        {
            const e2e::FiberFrameScan sc = e2e::scan_fiber_frames(dir);
            e2e::check(sc.frames == 3, "scan: кадров " + std::to_string(sc.frames));
            e2e::check(sc.health && sc.nonfinite == 0 && sc.first_problem.empty(),
                       "scan: health = true, нет не-конечных значений" +
                           (sc.first_problem.empty() ? std::string() : " (" + sc.first_problem + ")"));
            bool counts = sc.n_particles.size() == 3 && sc.n_bonds.size() == 3 && sc.time.size() == 3;
            if (counts)
            {
                counts = sc.n_particles[0] == nbeads && sc.n_particles[2] == nbeads && sc.n_bonds[0] == wbonds[0] &&
                         sc.n_bonds[1] == wbonds[1] && sc.n_bonds[2] == wbonds[2] && close_rel(sc.time[1], times[1]);
            }
            e2e::check(counts, "scan: число бусин / связей и TIME по кадрам");
        }

        std::printf("\n[5] Испорченные кадры\n");
        {
            const std::string good = slurp(files[0]);

            // 5a: "nan" вместо первой координаты первой частицы
            const fs::path bad = dir / "bad_nan";
            fs::create_directories(bad);
            std::string text = good;
            const std::size_t pos = text.find("PARTICLES\n");
            bool patched = false;
            if (pos != std::string::npos)
            {
                const std::size_t line0 = pos + std::string("PARTICLES\n").size();
                std::size_t p = text.find(' ', line0);          // после id
                p = text.find(' ', p + 1);                      // после object_id
                p = text.find(' ', p + 1);                      // после material_id -> перед x
                const std::size_t xe = text.find(' ', p + 1);   // конец x
                if (p != std::string::npos && xe != std::string::npos)
                {
                    text.replace(p + 1, xe - p - 1, "nan");
                    patched = true;
                }
            }
            e2e::check(patched, "испорченный кадр подготовлен (x первой частицы -> nan)");
            spit(bad / "fiber_output_00001.txt", text);
            const e2e::FiberFrameScan sb = e2e::scan_fiber_frames(bad);
            e2e::check(sb.frames == 1 && !sb.health && sb.nonfinite > 0 && !sb.first_problem.empty(),
                       "scan: nan обнаружен, nonfinite = " + std::to_string(sb.nonfinite) + ", " + sb.first_problem);
            e2e::FiberFrame fr;
            std::string err;
            e2e::check(e2e::read_fiber_frame(bad / "fiber_output_00001.txt", fr, &err) && std::isnan(fr.x[0][0]),
                       "nan не ломает разбор формата, читается как NaN");

            // 5b: обрезанный файл (меньше строк BONDS, чем в заголовке)
            const fs::path cut = dir / "bad_cut";
            fs::create_directories(cut);
            const std::size_t last_nl = good.find_last_of('\n', good.size() - 2);
            spit(cut / "fiber_output_00001.txt", good.substr(0, last_nl + 1));
            const e2e::FiberFrameScan sc = e2e::scan_fiber_frames(cut);
            err.clear();
            e2e::check(!e2e::read_fiber_frame(cut / "fiber_output_00001.txt", fr, &err) && !err.empty(),
                       "обрезанный файл: строгий разбор отказал (" + err + ")");
            e2e::check(!sc.health && !sc.first_problem.empty(), "scan: обрезанный кадр -> health = false");

            // 5c: лишняя строка данных сверх N_PARTICLES
            const fs::path extra = dir / "bad_extra";
            fs::create_directories(extra);
            std::string t3 = good;
            const std::size_t bp = t3.find("# i j l0");
            if (bp != std::string::npos) t3.insert(bp, "99 0 0 0 0 0 0 0 0\n");
            spit(extra / "fiber_output_00001.txt", t3);
            err.clear();
            e2e::check(!e2e::read_fiber_frame(extra / "fiber_output_00001.txt", fr, &err),
                       "лишняя строка частицы: счётчик не сошёлся (" + err + ")");

            // 5d: несуществующий файл и пустой каталог
            err.clear();
            e2e::check(!e2e::read_fiber_frame(dir / "fiber_output_99999.txt", fr, &err) && !err.empty(),
                       "несуществующий файл: ошибка чтения");
            const fs::path empty = dir / "empty_dir";
            fs::create_directories(empty);
            const e2e::FiberFrameScan se = e2e::scan_fiber_frames(empty);
            e2e::check(se.frames == 0 && e2e::list_fiber_frames(empty).empty(), "пустой каталог: 0 кадров");
        }
    } // ScopedWorkDir: временный каталог удалён, cwd восстановлен

    std::printf("\n[6] Побочные эффекты\n");
    e2e::check(fs::current_path() == cwd, "рабочий каталог восстановлен");
    e2e::check(e2e::list_fiber_frames(cwd).size() == fiber_files_before,
               "fiber_output_*.txt в рабочем каталоге не появились");
    e2e::check_no_side_effects(cwd, before);

    std::printf("\n=== e2e_fiber_frames: %s ===\n", e2e::failures == 0 ? "ПРОЙДЕН" : "ПРОВАЛЕН");
    return e2e::failures == 0 ? 0 : 1;
}
