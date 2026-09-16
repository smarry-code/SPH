// linked_list.cpp
// Определения функций поиска соседей на сетке (link_list) и замера времени
// выполнения программы (time_elapsed). Объявления — в linked_list.hpp.
#include "linked_list.hpp"
#include <cmath>
#include <iostream>
#include <omp.h>

// SetConsoleOutputCP (UTF-8 в Windows) объявлен в <windows.h>; в режиме
// разделения .hpp/.cpp каждый модуль — самостоятельная TU, поэтому windows.h
// подключается явно там, где используется консольный вывод.
#ifdef _WIN32
#include <windows.h>
#endif

namespace
{
// Статическая переменная для хранения времени начала выполнения программы
const auto program_start = std::chrono::high_resolution_clock::now();
}

void time_elapsed(double& elapsed)
{
    // Получаем текущее время
    auto now = std::chrono::high_resolution_clock::now();
    // Вычисляем разницу между текущим временем и временем старта программы
    std::chrono::duration<double> duration = now - program_start;
    // Сохраняем результат в переданную переменную
    elapsed = duration.count();
}

void link_list(const SimulationConfig& cfg, ParticleSystem& sys, size_t itimestep, size_t ntotal, size_t& niac, std::vector<size_t>& pair_i,
               std::vector<size_t>& pair_j, std::vector<double>& w, std::vector<std::vector<double>>& dwdx)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001); // Активация UTF-8 в Windows
#endif
    // === Замеры производительности ===
    double t_build = 0.0, t_pair = 0.0;
    double t_start, t_end;
    time_elapsed(t_start);

    int scale_k = 0;
    if (cfg.skf == 1) scale_k = 2;  // Кубический сплайн
    else if (cfg.skf == 2) scale_k = 3;  // Гауссово ядро
    else if (cfg.skf == 3) scale_k = 3;  // Квинтическое ядро

    // Этап 2: построение FlatGrid (замена CellMap)
    // Последовательная версия — оптимальна при N <= ~10000
    FlatGrid grid;
    grid.build(sys, ntotal, scale_k);
    double dx = grid.dx;

    niac = 0;
    time_elapsed(t_end);
    t_build = t_end - t_start;
    time_elapsed(t_start);

    // === Парный цикл (Этап 3: единый цикл + thread-local буфер dxiac) ===
    #pragma omp parallel
    {
        // thread-local переиспользуемый буфер dxiac.
        // Объявляется ОДИН раз на поток, переиспользуется без пересоздания.
        // std::array<double, dim> — стековый буфер точного размера, без аллокации кучи.
        std::array<double, dim> dxiac{};

        #pragma omp for schedule(auto)
        for (size_t i = 0; i < ntotal; ++i)
        {
            // Индекс ячейки частицы i
            int ci = grid.cell_index(
                sys.x[0][i],
                (dim >= 2) ? sys.x[1][i] : 0.0,
                (dim >= 3) ? sys.x[2][i] : 0.0);

            // Декомпозиция cell_idx → (ix, iy, iz)
            int ix = ci % grid.nx;
            int iy = (grid.nx > 0) ? (ci / grid.nx) % grid.ny : 0;
            int iz = (grid.nx * grid.ny > 0) ? ci / (grid.nx * grid.ny) : 0;

            int fac = static_cast<int>(std::round(scale_k * sys.hsml[i] / dx));

            // Сохранение семантики fac < 1: только собственная ячейка
            int rx = (fac < 1) ? 0 : fac;
            int ry = (dim >= 2) ? ((fac < 1) ? 0 : fac) : 0;
            int rz = (dim >= 3) ? ((fac < 1) ? 0 : fac) : 0;

            for (int kx = -rx; kx <= rx; ++kx)
            {
                int nx_cell = ix + kx;
                if (nx_cell < 0 || nx_cell >= grid.nx) continue;

                for (int ky = -ry; ky <= ry; ++ky)
                {
                    int ny_cell = iy + ky;
                    if (ny_cell < 0 || ny_cell >= grid.ny) continue;

                    for (int kz = -rz; kz <= rz; ++kz)
                    {
                        int nz_cell = iz + kz;
                        if (nz_cell < 0 || nz_cell >= grid.nz) continue;

                        int nc = nx_cell + grid.nx * (ny_cell + grid.ny * nz_cell);
                        auto [begin, end] = grid.cell_range(nc);

                        for (size_t idx = begin; idx < end; ++idx)
                        {
                            size_t j = grid.sorted_idx[idx];
                            if (j <= i) continue;

                            // Перезаписываем thread-local буфер (без аллокации)
                            double driac = 0.0;
                            for (size_t d = 0; d < dim; ++d)
                            {
                                dxiac[d] = sys.x[d][i] - sys.x[d][j];
                                driac += dxiac[d] * dxiac[d];
                            }
                            double mhsml = 0.5 * (sys.hsml[i] + sys.hsml[j]);
                            double cutoff = scale_k * mhsml;

                            // Проверка на взаимодействие
                            if (driac < cutoff * cutoff)
                            {
                                double r = std::sqrt(driac);
                                size_t pidx;
                                #pragma omp atomic capture
                                {
                                    pidx = niac;
                                    niac++;
                                }
                                if (pidx < cfg.max_interaction)
                                {
                                    double tdwdx[dim];
                                    double w_buf = 0.0;
                                    kernel(cfg.skf, r, dxiac, mhsml, w_buf, tdwdx);
                                    pair_i[pidx] = i;  // 0-базовая индексация
                                    pair_j[pidx] = j;
                                    w[pidx] = w_buf;
                                    for (size_t d = 0; d < dim; ++d)
                                    {
                                        dwdx[d][pidx] = tdwdx[d];
                                    }
                                }
                                else
                                {
                                    std::cout << "Превышен max_interaction, niac = " << niac << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    time_elapsed(t_end);
    t_pair = t_end - t_start;
    /*std::wcout << "[link_list] Построение сетки (FlatGrid): " << t_build << " с" << std::endl;
    std::wcout << "[link_list] Парный цикл: " << t_pair << " с" << std::endl;
    std::wcout << "[link_list] Итого: " << t_build + t_pair << " с" << std::endl;
    std::wcout << "[link_list] niac = " << niac << std::endl;
    std::wcout << "[link_list] n_cells = " << grid.n_cells() << " (" << grid.nx << "x" << grid.ny << "x" << grid.nz << ")" << std::endl;*/
}
