#pragma once
// diagnostics.hpp
// Единая валидация состояния частиц
// Назначение: единая точка валидации, вызываемая оркестратором single_step
// после ключевых этапов шага. Read-only + вывод в stderr: НЕ изменяет поля sys
// (dvxdt/dedt/rho/...) => численно нейтральна (побочные эффекты — только
// диагностика/лог). Существующие inline-проверки в силовых модулях сохранены
// как defence-in-depth (локальный контекст + управляемый отказ throw).
//
// Связь с управляемым отказом: validate возвращает bool; при !ok
// возможны счётчик аварий, статус шага, адаптивное уменьшение dt (полная
// схема — отдельная задача Приоритета 1).
#include "particle_system.hpp"
#include "rhs/cpu_executor.hpp"
#include <cmath>
#include <iostream>
#include <string>

struct Diagnostics
{
    // Проверяет конечность ключевых термодинамических полей для частиц [0, ntotal).
    // Возвращает true, если все значения конечны; иначе печатывает каждое
    // NaN/Inf в stderr и возвращает false.
    static bool validate(const ParticleSystem& sys, size_t ntotal, const std::string& stage)
    {
        bool ok = true;
        for (size_t i = 0; i < ntotal; ++i)
        {
            auto bad = [&](const char* f, double v)
            {
                if (!std::isfinite(v))
                {
                    std::cerr << "[NaN/Inf " << stage << "] " << f << "[" << i << "]=" << v << "\n";
                    ok = false;
                }
            };
            bad("p",    sys.p[i]);
            bad("rho",  sys.rho[i]);
            bad("c",    sys.c[i]);
            bad("u",    sys.u[i]);
            bad("hsml", sys.hsml[i]);
        }
        return ok;
    }

    // То же с исполнителем: параллельная проверка "есть ли не конечное значение" (булево ИЛИ -- точно при любом
    // порядке); сообщения печатает прежний последовательный проход, только если нарушение найдено.
    static bool validate(const rhs::CpuExecutor& exec, const ParticleSystem& sys, size_t ntotal,
                         const std::string& stage)
    {
        const int any_bad = exec.parallel_reduce(
            ntotal, 0,
            [&](size_t i) {
                return (!std::isfinite(sys.p[i]) || !std::isfinite(sys.rho[i]) || !std::isfinite(sys.c[i]) ||
                        !std::isfinite(sys.u[i]) || !std::isfinite(sys.hsml[i]))
                           ? 1
                           : 0;
            },
            [](int a, int b) { return a | b; });
        if (!any_bad) return true;
        return validate(sys, ntotal, stage);
    }
};
