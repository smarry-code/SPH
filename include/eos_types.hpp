#pragma once
// eos_types.hpp
// Чистые структуры параметров новых УРС (Linear, Shock). Без зависимости от
// Materials: заголовок может включаться и из eos.hpp, и из materials.hpp
// без цикла (ADR-005: УРС — чистые функции от параметров).

// Autodyn Linear: p = K·μ, μ = ρ/ρ0 − 1 (от энергии не зависит).
struct LinearEos
{
    double K = 0.0; // объёмный модуль, Па
};

// Autodyn Shock: Ми-Грюнайзен по адиабате Us = c0 + S1·up.
struct ShockEos
{
    double gamma0 = 0.0; // параметр Грюнайзена Γ0
    double c0 = 0.0;     // объёмная скорость звука, м/с
    double s1 = 0.0;     // наклон ударной адиабаты S1
};
