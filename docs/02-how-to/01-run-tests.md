# Запуск тестов и бенчмарков

В проекте настроены тесты через CTest (см. [`CMakeLists.txt`](../../../CMakeLists.txt:38)).

## Запуск всех тестов

```bash
ctest --test-dir build --output-on-failure
```

## Состав тестов

| Тест | Исполняемый файл | Что проверяет |
|---|---|---|
| `flat_grid_test` | `test_flat_grid` | Сеточный поиск соседей |
| `pair_graph_test` | `test_pair_graph` | Графовую раскраску пар (с OpenMP) |
| `particle_system_test` | `test_particle_system` | Тождество моста AoS↔SoA (точность 1e-12) |

## Бенчмарк

`bench_color` сравнивает последовательный (atomic) и параллельный (coloring)
пути обхода пар. Собирается отдельно и запускается вручную:

```bash
cmake --build build --target bench_color
./build/bench_color
```
