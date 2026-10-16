// fiber_output.cpp
// Вывод состояния нитей в файл. Объявление — в fiber/fiber_output.hpp.
//
// ФОРМАТ ФАЙЛА (fiber_output_<itimestep>.txt, текстовый):
//   строка 1: "TIME <t> STEP <itimestep> N_PARTICLES <n> N_BONDS <nb>"
//   затем блок PARTICLES: заголовок-комментарий, слово "PARTICLES", затем n
//     строк "<id> <object_id> <material_id> <x> <y> <z> <vx> <vy> <vz>"
//     (id — индекс бусины внутри FiberState, 0-based);
//   затем блок BONDS: заголовок-комментарий, слово "BONDS", затем nb строк
//     "<i> <j> <l0> <l> <l_break> <elongation>", где l — текущая длина связи
//     (i,j), elongation = (l - l0) / l0 (fiber/SPH.md, Fiber-Bond: l > l0 —
//     растяжение). l выводится тоже — оба варианта из пункта 4 плана.
//
//   РАСШИРЕНИЕ (неупругие материалы, fiber/fiber-plasticity-design.md, п. 10.3):
//   если среди материалов модели есть хоть один с неупругой связью
//   (FiberMaterial::inelastic()), строки связей и их комментарий-заголовок
//   получают четыре колонки В КОНЦЕ:
//     "# i j l0 l l_break elongation eps_p sigma w_p w_v"
//     eps_p — пластическая деформация связи (безразмерная);
//     sigma — напряжение связи [Па], ровно то, что считает RHS-ядро:
//       fiber::bond_stress(ε = l/l0 - 1, ε̇ = ((v_j - v_i)·n)/l0) от СОХРАНЁННОГО
//       состояния связи; скорости — финальные скорости бусин на момент вывода
//       (RHS берёт скорости предиктора, поэтому вязкий член может отличаться
//       на O(dt)). Сила в связи = A·sigma, A = pi d^2/4;
//     w_p, w_v — накопленная пластическая работа и вязкая диссипация связи
//       [Дж] (= BondState::w·A·l0).
//   Для упругих связей смешанной модели эти колонки: eps_p = 0,
//   sigma = E·max(ε, 0), w_p = w_v = 0. Если неупругих материалов нет (в том
//   числе при вызове без материалов), формат ровно прежний, 6 колонок --
//   читатели с жёсткой проверкой числа полей (e2e_fiber_frames.hpp) и
//   позиционные Python-читатели (tools/, experiments/) не ломаются.
#include "fiber/fiber_output.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "interactions/fiber_bond.hpp"

namespace fiber
{

void fiber_output(const Model& model, std::size_t itimestep, double time)
{
    fiber_output(model.fiber, model.topology, model.materials, itimestep, time);
}

void fiber_output(const FiberState& fiber, const FiberTopology& topology, std::size_t itimestep, double time)
{
    fiber_output(fiber, topology, std::vector<FiberMaterial>{}, itimestep, time);
}

void fiber_output(const FiberState& fiber, const FiberTopology& topology, const std::vector<FiberMaterial>& materials,
                  std::size_t itimestep, double time)
{
    // Формат x[dim]/vx[dim] осмыслен только при dim == 3 (см.
    // check_fiber_requires_dim3 в fiber_system.hpp); проверяем и здесь,
    // отдельно от load_fiber_file, чтобы вызов вывода в неверной сборке не
    // писал молча бессмысленный файл.
    check_fiber_requires_dim3();

    std::ostringstream filename;
    filename << "fiber_output_" << std::setw(5) << std::setfill('0') << itimestep << ".txt";

    std::ofstream file(filename.str());
    if (!file.is_open())
    {
        throw std::runtime_error("fiber::fiber_output: не удалось открыть файл для вывода данных нитей");
    }

    const std::size_t n = fiber.size();
    const std::size_t nb = topology.bonds.size();

    file << "TIME " << time << " STEP " << itimestep
         << " N_PARTICLES " << n << " N_BONDS " << nb << "\n";

    file << "# id object_id material_id x y z vx vy vz\n";
    file << "PARTICLES\n";
    for (std::size_t i = 0; i < n; ++i)
    {
        file << i << " " << fiber.object_id[i] << " " << fiber.material_id[i];
        for (std::size_t d = 0; d < dim; ++d) file << " " << fiber.x[d][i];
        for (std::size_t d = 0; d < dim; ++d) file << " " << fiber.vx[d][i];
        file << "\n";
    }

    bool extended = false;
    for (const auto& m : materials) extended = extended || m.inelastic();

    file << (extended ? "# i j l0 l l_break elongation eps_p sigma w_p w_v\n" : "# i j l0 l l_break elongation\n");
    file << "BONDS\n";
    for (std::size_t k = 0; k < nb; ++k)
    {
        const FiberBond& b = topology.bonds[k];
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dc = fiber.x[d][b.j] - fiber.x[d][b.i];
            l2 += dc * dc;
        }
        const double l = std::sqrt(l2);
        const double elongation = (b.l0 > 0.0) ? (l - b.l0) / b.l0 : 0.0;
        file << b.i << " " << b.j << " " << b.l0 << " " << l << " " << b.l_break << " " << elongation;
        if (extended)
        {
            const FiberMaterial& mat = materials[static_cast<std::size_t>(fiber.material_id[b.i])];
            const double A = interactions::fiber_cross_section_area(mat);
            const double eps = (b.l0 > 0.0) ? l / b.l0 - 1.0 : 0.0;
            double eps_p = 0.0, sigma = 0.0, wp = 0.0, wv = 0.0;
            if (mat.inelastic())
            {
                double eps_dot = 0.0;
                const ViscoParams vp = interactions::fiber_bond_visco(mat, A, b.l0, fiber.mass[b.i], fiber.mass[b.j]);
                if (vp.eta != 0.0 && l > 0.0)
                {
                    double vrel = 0.0;
                    for (std::size_t d = 0; d < dim; ++d)
                        vrel += (fiber.vx[d][b.j] - fiber.vx[d][b.i]) * (fiber.x[d][b.j] - fiber.x[d][b.i]);
                    eps_dot = (vrel / l) / b.l0;
                }
                const BondLaw law(mat.E, mat.plastic, mat.rate, vp, mat.failure);
                eps_p = b.st.eps_p;
                sigma = bond_stress(law, b.st, eps, eps_dot);
                wp = b.st.w_p * A * b.l0;
                wv = b.st.w_v * A * b.l0;
            }
            else
            {
                sigma = mat.E * std::max(eps, 0.0);
            }
            // Новые колонки -- с повышенной точностью (17 значащих цифр:
            // eps_p ~ 1e-5 и разность sigma ~ E*eps при 6 цифрах теряются);
            // прежние 6 колонок и точность потока не трогаем.
            const std::streamsize old_prec = file.precision(17);
            file << " " << eps_p << " " << sigma << " " << wp << " " << wv;
            file.precision(old_prec);
        }
        file << "\n";
    }

    file.close();
}

} // namespace fiber
