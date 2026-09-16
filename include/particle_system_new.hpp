#pragma once

#include "param.hpp"

#include <vector>

typedef float fp;

template <typename T, std::size_t D>
using vecs = std::array<T, D>;

struct ParticleSystem
{
    HydrodynamicParticles hydrodynamic_particles;
    FiberParticles fiber_particles;
    RigidbodyParticles rigidbody_particles;

    //TODO: methods
    //Методы служебные, на данный момент не важно
};

struct Dynamics
{
    vecs<fp, dim> x;
    vecs<fp, dim> vx;
};


struct Thermodynamics
{
    std::vector<fp> mass, rho, rho_0, p, u, c, s, e, hsml, dp, e_pl, D;
};


struct HydrodynamicParticles
{
    Dynamics dynamics;
    Thermodynamics thermodynamics;
};


struct FiberParticles
{
    Dynamics dynamics;
    Thermodynamics thermodynamics;

    //TODO: fiber topology
    std::vector<size_t> fiber_offsets;  // границы нитей

    std::vector<fp> rest_length;      // длина для каждого сегмента

    // std::vector<fp> young_modulus; в thermodynamics   // модуль Юнга
    std::vector<fp> diameter;         // диаметр
};


struct RigidbodyParticles
{
    Dynamics dynamics;
    Thermodynamics thermodynamics;

    // --- тензор напряжений (девиатор) ---
    std::vector<double> sxx, syy, szz, sxy, syz, sxz;

    // --- тензор скоростей деформаций ---
    std::vector<double> exx, eyy, ezz, exy, exz, eyz;

    // --- тензор скоростей поворота ---
    std::vector<double> rxy, rxz, ryz;
};


struct LocalPrimitive
{
    std::vector<size_t> get_particles();

    // TODO: return type
    void compute_interaction(LocalPrimitive& other);

};


struct HydrodynamicLocalPrimitive : public LocalPrimitive
{
    HydrodynamicLocalPrimitive(size_t particle_idx);
    // LocalPrimitive is single particle
};

struct FiberLocalPrimitive : public LocalPrimitive
{
    FiberLocalPrimitive(size_t particle_idx);
    // LocalPrimitive is [particle_idx - 1, particle_idx, particle_idx + 1]
};

struct RigidbodyLocalPrimitive : public LocalPrimitive
{
    RigidbodyLocalPrimitive(size_t particle_idx);
    //TODO: разобраться
};

void compute_interaction(HydrodynamicLocalPrimitive& prim_i, HydrodynamicLocalPrimitive& prim_j);
void compute_interaction(HydrodynamicLocalPrimitive& prim_i, FiberLocalPrimitive& prim_j);
void compute_interaction(HydrodynamicLocalPrimitive& prim_i, RigidbodyLocalPrimitive& prim_j);

void compute_interaction(FiberLocalPrimitive& prim_i, HydrodynamicLocalPrimitive& prim_j);
void compute_interaction(FiberLocalPrimitive& prim_i, FiberLocalPrimitive& prim_j);
void compute_interaction(FiberLocalPrimitive& prim_i, RigidbodyLocalPrimitive& prim_j);

void compute_interaction(RigidbodyLocalPrimitive& prim_i, HydrodynamicLocalPrimitive& prim_j);
void compute_interaction(RigidbodyLocalPrimitive& prim_i, FiberLocalPrimitive& prim_j);
void compute_interaction(RigidbodyLocalPrimitive& prim_i, RigidbodyLocalPrimitive& prim_j);

