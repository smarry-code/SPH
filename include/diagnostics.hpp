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
            auto bad = [&](const char* f, real_t v)
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
};
