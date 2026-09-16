// virt_part.cpp
// Определение функции генерации виртуальных частиц.
// Объявление — в virt_part.hpp.
#include "virt_part.hpp"

void virt_part(ParticleSystem& sys, size_t ntotal, size_t& nvirt)
{

    // Инициализация nvirt
    nvirt = 0;
    if (axis_sim == true && dim == 2)
    {
        for (size_t i = 0; i < ntotal; ++i)
        {
            if (sys.x[0][i] < 2 * sys.hsml[i])
            {
                sys.inzializ_2d_symm_index(nvirt + ntotal, i);
                nvirt++;
            }
        }
    }
    // Сохранение данных о виртуальных частицах
    /*if (nvirt > 0 && (itimestep % save_step == 0))
    {
        std::ofstream xv_file("../data/xy_vp.dat");
        std::ofstream state_file("../data/state_vp.dat");
        std::ofstream other_file("../data/other_vp.dat");

        if (!xv_file.is_open() || !state_file.is_open() || !other_file.is_open())
        {
            throw std::runtime_error("Не удалось создать файлы для сохранения виртуальных частиц");
        }

        xv_file << nvirt << "\n";

        // Запись данных
        for (size_t i = ntotal; i < ntotal + nvirt; ++i)
        {
            // Форматированный вывод координат и скоростей
            xv_file << i + 1; // Индекс в формате Фортрана (начинается с 1)
            for (size_t d = 0; d < dim; ++d)
            {
                xv_file << " " << std::scientific << sys.x[d][i];
            }
            for (size_t d = 0; d < dim; ++d)
            {
                xv_file << " " << std::scientific << sys.vx[d][i];
            }
            xv_file << "\n";

            // Форматированный вывод состояния
            state_file << i + 1 << " "
                       << std::scientific << sys.mass[i] << " "
                       << std::scientific << sys.rho[i] << " "
                       << std::scientific << sys.p[i] << " "
                       << std::scientific << sys.u[i] << "\n";

            // Форматированный вывод других параметров
            other_file << i + 1 << " "
                       << sys.itype[i] << " "
                       << std::scientific << sys.hsml[i] << "\n";
        }

        xv_file.close();
        state_file.close();
        other_file.close();
    }

    // Вывод статистики
    if (int_stat && itimestep % print_step == 0)
    {
        std::cout << " >> Статистика: виртуальные граничные частицы\n";
        std::cout << "    Число виртуальных частиц: " << nvirt << "\n";
    }*/

}
