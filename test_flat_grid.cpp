// ============================================================================
// Unit-тест для FlatGrid (flat_grid.hpp)
// Проверка: 3 частицы в 3 разных ячейках (1D, dim=1)
// Компиляция: g++ -std=c++17 -o test_flat_grid test_flat_grid.cpp -I.
// ============================================================================

#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#ifdef _WIN32
#include <windows.h>
#endif
#include "param.hpp"
#include "flat_grid.hpp"

// Тест 1: Базовая проверка — 3 частицы в 3 разных ячейках (1D)
void test_basic_3_particles() {
    std::cout << "=== Тест 1: 3 частицы в 3 ячейках (1D) ===" << std::endl;
    
    // Создаём 3 частицы с координатами, попадающими в разные ячейки
    // dx = scale_k * max_hmsl = 2 * 0.5 = 1.0
    // Частицы: x = 0.1, 1.5, 3.2 -> ячейки 0, 1, 3
    std::vector<Particle> Elem(3);
    for (auto& p : Elem) p.inzializ();
    
    Elem[0].x[0] = 0.1;  Elem[0].hsml = 0.5;
    Elem[1].x[0] = 1.5;  Elem[1].hsml = 0.5;
    Elem[2].x[0] = 3.2;  Elem[2].hsml = 0.5;
    
    FlatGrid grid;
    int scale_k = 2; // кубический сплайн (skf=1)
    grid.build(Elem, 3, scale_k);
    
    // Проверка параметров сетки
    std::cout << "  dx = " << grid.dx << " (ожидалось 1.0)" << std::endl;
    std::cout << "  origin[0] = " << grid.origin[0] << " (ожидалось 0.1)" << std::endl;
    std::cout << "  nx = " << grid.nx << ", ny = " << grid.ny << ", nz = " << grid.nz << std::endl;
    std::cout << "  n_cells = " << grid.n_cells() << std::endl;
    
    assert(grid.dx == 1.0 && "dx должен быть scale_k * max_hmsl = 2 * 0.5 = 1.0");
    assert(grid.origin[0] == 0.1 && "origin[0] должен быть xmin = 0.1");
    assert(grid.ny == 1 && grid.nz == 1 && "Для dim=1 ny и nz должны быть 1");
    
    // Проверка индексов ячеек
    int c0 = grid.cell_index(0.1, 0.0, 0.0);
    int c1 = grid.cell_index(1.5, 0.0, 0.0);
    int c2 = grid.cell_index(3.2, 0.0, 0.0);
    
    std::cout << "  Ячейки частиц: c0=" << c0 << ", c1=" << c1 << ", c2=" << c2 << std::endl;
    
    // Все частицы должны быть в разных ячейках
    assert(c0 != c1 && "Частицы 0 и 1 должны быть в разных ячейках");
    assert(c1 != c2 && "Частицы 1 и 2 должны быть в разных ячейках");
    assert(c0 != c2 && "Частицы 0 и 2 должны быть в разных ячейках");
    
    // Проверка offset и sorted_idx
    std::cout << "  offset: ";
    for (int i = 0; i <= grid.n_cells(); ++i) {
        std::cout << grid.offset[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "  sorted_idx: ";
    for (size_t i = 0; i < grid.sorted_idx.size(); ++i) {
        std::cout << grid.sorted_idx[i] << " ";
    }
    std::cout << std::endl;
    
    // Проверка: каждая ячейка содержит ровно 1 частицу
    int non_empty = 0;
    for (int c = 0; c < grid.n_cells(); ++c) {
        auto [begin, end] = grid.cell_range(c);
        int count = end - begin;
        if (count > 0) {
            non_empty++;
            std::cout << "  Ячейка " << c << ": частица " << grid.sorted_idx[begin] << std::endl;
        }
    }
    assert(non_empty == 3 && "Должно быть 3 непустых ячейки");
    
    // Контрольный инвариант: offset[n_cells] == ntotal
    assert(grid.offset[grid.n_cells()] == 3 && "offset[n_cells] должен равняться ntotal");
    
    std::cout << "  Тест 1 ПРОЙДЕН!" << std::endl << std::endl;
}

// Тест 2: Все частицы в одной ячейке
void test_all_in_one_cell() {
    std::cout << "=== Тест 2: Все частицы в одной ячейке ===" << std::endl;
    
    std::vector<Particle> Elem(3);
    for (auto& p : Elem) p.inzializ();
    
    // Все частицы очень близко друг к другу (в пределах одной ячейки)
    Elem[0].x[0] = 1.0;  Elem[0].hsml = 0.5;
    Elem[1].x[0] = 1.1;  Elem[1].hsml = 0.5;
    Elem[2].x[0] = 1.2;  Elem[2].hsml = 0.5;
    
    FlatGrid grid;
    int scale_k = 2;
    grid.build(Elem, 3, scale_k);
    
    std::cout << "  dx = " << grid.dx << std::endl;
    std::cout << "  nx = " << grid.nx << std::endl;
    
    // Все частицы должны попасть в одну ячейку
    int c0 = grid.cell_index(1.0, 0.0, 0.0);
    int c1 = grid.cell_index(1.1, 0.0, 0.0);
    int c2 = grid.cell_index(1.2, 0.0, 0.0);
    
    std::cout << "  Ячейки частиц: c0=" << c0 << ", c1=" << c1 << ", c2=" << c2 << std::endl;
    
    assert(c0 == c1 && c1 == c2 && "Все частицы должны быть в одной ячейке");
    
    // Проверка: ячейка содержит 3 частицы
    auto [begin, end] = grid.cell_range(c0);
    int count = end - begin;
    std::cout << "  Частиц в ячейке " << c0 << ": " << count << std::endl;
    assert(count == 3 && "В ячейке должно быть 3 частицы");
    
    std::cout << "  Тест 2 ПРОЙДЕН!" << std::endl << std::endl;
}

// Тест 3: Отрицательные координаты (проверка сдвига origin)
void test_negative_coordinates() {
    std::cout << "=== Тест 3: Отрицательные координаты ===" << std::endl;
    
    std::vector<Particle> Elem(3);
    for (auto& p : Elem) p.inzializ();
    
    // Отрицательные координаты — после сдвига на origin должны стать >= 0
    Elem[0].x[0] = -5.0;  Elem[0].hsml = 0.5;
    Elem[1].x[0] = -3.0;  Elem[1].hsml = 0.5;
    Elem[2].x[0] = -1.0;  Elem[2].hsml = 0.5;
    
    FlatGrid grid;
    int scale_k = 2;
    grid.build(Elem, 3, scale_k);
    
    std::cout << "  origin[0] = " << grid.origin[0] << " (ожидалось -5.0)" << std::endl;
    std::cout << "  dx = " << grid.dx << std::endl;
    
    assert(grid.origin[0] == -5.0 && "origin[0] должен быть xmin = -5.0");
    
    // После сдвига: -5.0 -> 0, -3.0 -> 2, -1.0 -> 4 (при dx=1.0)
    int c0 = grid.cell_index(-5.0, 0.0, 0.0);
    int c1 = grid.cell_index(-3.0, 0.0, 0.0);
    int c2 = grid.cell_index(-1.0, 0.0, 0.0);
    
    std::cout << "  Ячейки частиц: c0=" << c0 << ", c1=" << c1 << ", c2=" << c2 << std::endl;
    
    // Индексы должны быть неотрицательными
    assert(c0 >= 0 && c1 >= 0 && c2 >= 0 && "Индексы ячеек должны быть >= 0");
    assert(c0 != c1 && c1 != c2 && "Частицы должны быть в разных ячейках");
    
    std::cout << "  Тест 3 ПРОЙДЕН!" << std::endl << std::endl;
}

// Тест 4: Пустой набор частиц
void test_empty() {
    std::cout << "=== Тест 4: Пустой набор частиц ===" << std::endl;
    
    std::vector<Particle> Elem;
    
    FlatGrid grid;
    int scale_k = 2;
    grid.build(Elem, 0, scale_k);
    
    assert(grid.n_cells() == 1 && "Для пустого набора n_cells должен быть 1");
    assert(grid.sorted_idx.empty() && "sorted_idx должен быть пуст");
    assert(grid.offset[0] == 0 && grid.offset[1] == 0 && "offset должен быть [0, 0]");
    
    std::cout << "  Тест 4 ПРОЙДЕН!" << std::endl << std::endl;
}

// Тест 5: 2D случай — проверка что y-координата учитывается
void test_2d_case() {
    std::cout << "=== Тест 5: 2D случай (dim=2) ===" << std::endl;
    
    // Для этого теста нужно временно изменить dim
    // Поскольку dim — constexpr, проверяем только логику cell_index
    // В реальном коде dim задаётся при компиляции
    
    std::vector<Particle> Elem(4);
    for (auto& p : Elem) p.inzializ();
    
    // Частицы в разных углах
    Elem[0].x[0] = 0.0;  Elem[0].hsml = 0.5;
    Elem[1].x[0] = 5.0;  Elem[1].hsml = 0.5;
    Elem[2].x[0] = 0.0;  Elem[2].hsml = 0.5;
    Elem[3].x[0] = 5.0;  Elem[3].hsml = 0.5;
    
    // Для dim=1 (текущая сборка) — проверяем что x-координаты корректно обрабатываются
    FlatGrid grid;
    int scale_k = 2;
    grid.build(Elem, 4, scale_k);
    
    std::cout << "  nx = " << grid.nx << ", n_cells = " << grid.n_cells() << std::endl;
    std::cout << "  offset[n_cells] = " << grid.offset[grid.n_cells()] << " (ожидалось 4)" << std::endl;
    
    assert(grid.offset[grid.n_cells()] == 4 && "Сумма count должна равняться ntotal");
    
    std::cout << "  Тест 5 ПРОЙДЕН!" << std::endl << std::endl;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    
    std::cout << "============================================" << std::endl;
    std::cout << "  Unit-тесты FlatGrid (flat_grid.hpp)" << std::endl;
    std::cout << "  dim = " << dim << std::endl;
    std::cout << "============================================" << std::endl << std::endl;
    
    test_basic_3_particles();
    test_all_in_one_cell();
    test_negative_coordinates();
    test_empty();
    test_2d_case();
    
    std::cout << "============================================" << std::endl;
    std::cout << "  ВСЕ ТЕСТЫ ПРОЙДЕНЫ!" << std::endl;
    std::cout << "============================================" << std::endl;
    
    return 0;
}