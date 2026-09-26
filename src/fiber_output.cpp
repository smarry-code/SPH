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
#include "fiber/fiber_output.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace fiber
{

void fiber_output(const Model& model, std::size_t itimestep, double time)
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

    const std::size_t n = model.fiber.size();
    const std::size_t nb = model.topology.bonds.size();

    file << "TIME " << time << " STEP " << itimestep
         << " N_PARTICLES " << n << " N_BONDS " << nb << "\n";

    file << "# id object_id material_id x y z vx vy vz\n";
    file << "PARTICLES\n";
    for (std::size_t i = 0; i < n; ++i)
    {
        file << i << " " << model.fiber.object_id[i] << " " << model.fiber.material_id[i];
        for (std::size_t d = 0; d < dim; ++d) file << " " << model.fiber.x[d][i];
        for (std::size_t d = 0; d < dim; ++d) file << " " << model.fiber.vx[d][i];
        file << "\n";
    }

    file << "# i j l0 l l_break elongation\n";
    file << "BONDS\n";
    for (std::size_t k = 0; k < nb; ++k)
    {
        const FiberBond& b = model.topology.bonds[k];
        double l2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double dc = static_cast<double>(model.fiber.x[d][b.j]) - static_cast<double>(model.fiber.x[d][b.i]);
            l2 += dc * dc;
        }
        const double l = std::sqrt(l2);
        const double l0 = static_cast<double>(b.l0); // вывод -- в double (precision.hpp)
        const double elongation = (l0 > 0.0) ? (l - l0) / l0 : 0.0;
        file << b.i << " " << b.j << " " << b.l0 << " " << l << " " << b.l_break << " " << elongation << "\n";
    }

    file.close();
}

} // namespace fiber
