// fiber_input.cpp
// Разбор файла ломаных нитей в fiber::Model (fiber/PLAN.md, этап Ф2, пункт 3).
// Объявление load_fiber_file — в include/fiber/fiber_input.hpp.
//
// ФОРМАТ ФАЙЛА (текстовый, свободный формат, разделители — пробелы):
//   - строки-комментарии и хвостовые комментарии начинаются с '#';
//   - пустые строки игнорируются;
//   - "MATERIAL <id> <E> <nu> <d> <rho> <eps_break> [mu_ff]"
//       определяет материал нити с индексом <id> (целое >= 0; индексы не
//       обязаны быть подряд, но должны быть объявлены ДО первой ссылки на
//       них из FIBER). E, nu, d, rho, eps_break — double (модуль Юнга,
//       коэффициент Пуассона, диаметр, плотность, критическое удлинение
//       разрыва; см. fiber/SPH.md, разделы Fiber-Bond и Fiber-Solid-Contact).
//       mu_ff — НЕОБЯЗАТЕЛЬНОЕ седьмое поле, коэффициент динамического
//       трения нить-нить (fiber/fiber_friction_implementation.md, §8);
//       если не задан, mu_ff = 0 (старое поведение без трения). Если задан,
//       должен быть числом >= 0, иначе — ошибка формата.
//   - "FIBER <material_id> <point_count>"
//       открывает одну ломаную (= одну нить). За строкой заголовка следуют
//       РОВНО <point_count> строк "x y z" — вершины ломаной в порядке вдоль
//       нити, в ненагруженной (недеформированной) конфигурации. <point_count>
//       >= 2. Каждая ломаная становится отдельной Fiber с собственным
//       object_id (0, 1, 2, ... в порядке следования блоков FIBER в файле).
//
// Пример:
//   # Нить из капрона, две прямые ломаные по 3 и 2 точки
//   MATERIAL 0 2.0e9 0.35 1.0e-4 1140.0 0.02
//   FIBER 0 3
//   0.0 0.0 0.0
//   0.01 0.0 0.0
//   0.02 0.0 0.0
//   FIBER 0 2
//   0.0 0.01 0.0
//   0.0 0.02 0.0
//
// Вычисления при разборе (fiber/SPH.md, раздел «Fiber»):
//   l0 связи (k, k+1) = расстояние между точками файла (как заданы, без
//   деформации);
//   l_break = (1 + eps_break) * l0 (материал связи; fiber/SPH.md,
//   Fiber-Bond, «Постановка проекта: только растяжение»);
//   масса бусины k: m_k = rho_f * A * (l_{k-1,k} + l_{k,k+1}) / 2,
//   A = pi*d^2/4, для концевой бусины берётся только один прилежащий
//   сегмент (второе слагаемое отсутствует, а не делится пополам дважды).
#include "fiber/fiber_input.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "param.hpp"

namespace fiber
{

namespace
{

struct RawFiber
{
    int material_id = -1;
    std::vector<std::array<double, 3>> points;
};

// Обрезает комментарий, начинающийся с '#'. Возвращает true, если после
// обрезки строка непустая (содержит не-пробельные символы).
bool strip_comment_nonblank(std::string& line)
{
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    const std::size_t b = line.find_first_not_of(" \t\r\n");
    return b != std::string::npos;
}

} // namespace

void load_fiber_file(const std::string& filename, Model& model)
{
    // Границу Р9 проверяем прежде всего: дальше код пишет ровно dim
    // компонент в FiberState::x[dim]/vx[dim] из троек (x,y,z) файла, и это
    // корректно только при dim == 3 (иначе часть заданных координат молча
    // терялась бы).
    check_fiber_requires_dim3();

    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("fiber::load_fiber_file: не удалось открыть файл '" + filename + "'");
    }

    std::vector<FiberMaterial> materials;
    std::vector<char> material_defined;
    std::vector<RawFiber> fibers;

    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line))
    {
        ++line_no;
        if (!strip_comment_nonblank(line)) continue;

        std::istringstream iss(line);
        std::string kw;
        iss >> kw;

        if (kw == "MATERIAL")
        {
            long long id = -1;
            double E = 0.0, nu = 0.0, d = 0.0, rho = 0.0, eps_break = 0.0;
            if (!(iss >> id >> E >> nu >> d >> rho >> eps_break) || id < 0)
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: неверный формат MATERIAL в строке " + std::to_string(line_no));
            }

            // mu_ff (коэффициент динамического трения нить-нить,
            // fiber/fiber_friction_implementation.md, §8) — НЕОБЯЗАТЕЛЬНОЕ
            // седьмое поле строки MATERIAL. Если после шести обязательных
            // полей в строке больше нет токена, оставляем mu_ff = 0.0
            // (старое поведение без трения, ничего не меняющее для уже
            // существующих файлов нитей). Если токен есть, он обязан
            // разбираться как double ЦЕЛИКОМ (иначе — ошибка формата, как и
            // для остальных полей MATERIAL выше) и быть неотрицательным.
            // Разбор — потоком, как у остальных полей (локаль потока
            // "classic"), а не std::stod, зависящим от глобальной C-локали.
            double mu_ff = 0.0;
            std::string mu_ff_tok;
            if (iss >> mu_ff_tok)
            {
                std::istringstream tok(mu_ff_tok);
                const bool ok = static_cast<bool>(tok >> mu_ff) && (tok >> std::ws).eof();
                if (!ok)
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: неверный формат MATERIAL (mu_ff) в строке " +
                        std::to_string(line_no));
                }
                if (mu_ff < 0.0)
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: mu_ff должен быть >= 0 в MATERIAL, строка " +
                        std::to_string(line_no));
                }
            }

            const std::size_t uid = static_cast<std::size_t>(id);
            if (uid >= materials.size())
            {
                materials.resize(uid + 1);
                material_defined.resize(uid + 1, 0);
            }
            materials[uid] = FiberMaterial{E, nu, d, rho, eps_break, mu_ff};
            material_defined[uid] = 1;
        }
        else if (kw == "FIBER")
        {
            long long material_id = -1;
            long long point_count = 0;
            if (!(iss >> material_id >> point_count) || material_id < 0 || point_count < 2)
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: неверный заголовок FIBER (материал/число точек) в строке " +
                    std::to_string(line_no));
            }
            const std::size_t umat = static_cast<std::size_t>(material_id);
            if (umat >= materials.size() || !material_defined[umat])
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: FIBER ссылается на неопределённый материал " +
                    std::to_string(material_id) + " в строке " + std::to_string(line_no));
            }

            RawFiber rf;
            rf.material_id = material_id;
            rf.points.reserve(static_cast<std::size_t>(point_count));
            for (long long p = 0; p < point_count; ++p)
            {
                if (!std::getline(in, line))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: неожиданный конец файла внутри FIBER (после строки " +
                        std::to_string(line_no) + ")");
                }
                ++line_no;
                strip_comment_nonblank(line);
                std::istringstream piss(line);
                std::array<double, 3> pt{};
                if (!(piss >> pt[0] >> pt[1] >> pt[2]))
                {
                    throw std::runtime_error(
                        "fiber::load_fiber_file: неверный формат точки в строке " + std::to_string(line_no));
                }
                rf.points.push_back(pt);
            }
            fibers.push_back(std::move(rf));
        }
        else
        {
            throw std::runtime_error(
                "fiber::load_fiber_file: неизвестное ключевое слово '" + kw + "' в строке " + std::to_string(line_no));
        }
    }

    if (fibers.empty())
    {
        throw std::runtime_error("fiber::load_fiber_file: файл '" + filename + "' не содержит ни одной нити (FIBER)");
    }

    std::size_t total_n = 0;
    for (const auto& f : fibers) total_n += f.points.size();

    model.materials = materials;
    model.fiber.resize(total_n);
    model.topology.resize(total_n);
    model.topology.bonds.clear();
    model.topology.bonds.reserve(total_n);
    ++model.topology.version;

    std::size_t base = 0;
    for (std::size_t fid = 0; fid < fibers.size(); ++fid)
    {
        const RawFiber& f = fibers[fid];
        const std::size_t n = f.points.size();
        const FiberMaterial& mat = model.materials[static_cast<std::size_t>(f.material_id)];
        const double A = pi * mat.d * mat.d / 4.0;

        std::vector<double> seg(n - 1);
        for (std::size_t k = 0; k + 1 < n; ++k)
        {
            double d2 = 0.0;
            for (std::size_t c = 0; c < 3; ++c)
            {
                const double dc = f.points[k + 1][c] - f.points[k][c];
                d2 += dc * dc;
            }
            seg[k] = std::sqrt(d2);
            if (!(seg[k] > 0.0))
            {
                throw std::runtime_error(
                    "fiber::load_fiber_file: нулевая или отрицательная длина сегмента в нити " +
                    std::to_string(fid) + " (точки " + std::to_string(k) + ", " + std::to_string(k + 1) + ")");
            }
        }

        for (std::size_t k = 0; k < n; ++k)
        {
            const std::size_t gi = base + k;
            for (std::size_t d = 0; d < dim; ++d)
            {
                model.fiber.x[d][gi] = f.points[k][d];
                model.fiber.vx[d][gi] = 0.0;
            }
            const double l_prev = (k > 0) ? seg[k - 1] : 0.0;
            const double l_next = (k + 1 < n) ? seg[k] : 0.0;
            model.fiber.mass[gi] = mat.rho * A * (l_prev + l_next) / 2.0;
            model.fiber.object_id[gi] = static_cast<int>(fid);
            model.fiber.material_id[gi] = f.material_id;

            model.topology.prev[gi] = (k > 0) ? static_cast<long long>(gi - 1) : -1;
            model.topology.next[gi] = (k + 1 < n) ? static_cast<long long>(gi + 1) : -1;
        }

        for (std::size_t k = 0; k + 1 < n; ++k)
        {
            FiberBond b;
            b.i = base + k;
            b.j = base + k + 1;
            b.l0 = seg[k];
            b.l_break = (1.0 + mat.eps_break) * seg[k];
            model.topology.bonds.push_back(b);
        }

        base += n;
    }
}

} // namespace fiber
