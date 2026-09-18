// buffers.hpp
// общий пул переиспользуемых рабочих буферов.
//
// Назначение: устранить аллокацию/обнуление массивов размера ntotal на каждом
// шаге интегрирования. Раньше функции силового расчёта (int_force и др.)
// создавали локальные std::vector<double>(ntotal, 0.0) при каждом вызове,
// что на длинных симуляциях приводило к миллионам malloc/free. WorkBuffers
// растёт только при увеличении ntotal, а перед использованием обнуляется в
// диапазоне [0, n) методом clear_tensor().
//
// Безопасность использования:
//   Экземпляр создаётся ОДИН раз на всю симуляцию (в time_integration) и
//   передаётся по ссылке в single_step -> int_force. Силовые функции
//   вызываются последовательно по одной за такт, а внутренний OpenMP-
//   параллелизм только читает тензорные буферы (запись выполняется в
//   последовательном блоке до входа в параллельную область). Поэтому
//   используется обычная (не thread_local) ссылка — разделяемое состояние
//   безопасно при отсутствии одновременных вызовов.
#pragma once
#include <vector>
#include <cstddef>
#include <algorithm>

struct WorkBuffers
{
    // Компоненты тензора скоростей сдвига (заполняются в int_force при visc=true,
    // читаются в параллельном цикле расчёта сил). Размер = ntotal.
    std::vector<double> txx;
    std::vector<double> tyy;
    std::vector<double> tzz;
    std::vector<double> txy;
    std::vector<double> txz;
    std::vector<double> tyz;

    // --- Контакт (ADR-007). Пишутся в ext_force, читаются в time_integration.
    // ext_force обходит пары последовательно (OpenMP там нет), гонок нет.
    double min_dt_contact    = 1e30;  // ограничитель шага; 1e30 => контакт не ограничивает (contact::kNoLimit)
    size_t contact_pair_count = 0;    // число пар с p > 0 за шаг
    double contact_max_delta  = 0.0;  // max(p) за шаг, м

    // Сброс контактной диагностики перед накоплением на текущем шаге.
    void reset_contact()
    {
        min_dt_contact = 1e30;
        contact_pair_count = 0;
        contact_max_delta = 0.0;
    }

    // Гарантирует вместимость n элементов во всех тензорных буферах.
    // Растёт только вверх (без усечения), поэтому на шагах после первого — no-op.
    void ensure(size_t n)
    {
        if (txx.size() < n)
        {
            txx.resize(n);
            tyy.resize(n);
            tzz.resize(n);
            txy.resize(n);
            txz.resize(n);
            tyz.resize(n);
        }
    }

    // Обнуляет первые n элементов тензорных буферов (приведение к исходному
    // состоянию перед накоплением/чтением на текущем шаге).
    void clear_tensor(size_t n)
    {
        std::fill_n(txx.data(), n, 0.0);
        std::fill_n(tyy.data(), n, 0.0);
        std::fill_n(tzz.data(), n, 0.0);
        std::fill_n(txy.data(), n, 0.0);
        std::fill_n(txz.data(), n, 0.0);
        std::fill_n(tyz.data(), n, 0.0);
    }
};
