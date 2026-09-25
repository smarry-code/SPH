// particle_system.hpp
// Structure-of-Arrays (SoA) представление данных частиц.
//
// Назначение: заменить std::vector<Particle> (AoS с двойной косвенностью:
// Elem[i].vx[d] = чтение указателя Particle + чтение указателя внутреннего
// буфера vx + промах кэша) на совокупность плоских непрерывных массивов.
// Доступ sys.rho[i] — ОДНО разыменование, данные смежных индексов лежат
// подряд => prefetch, векторная загрузка, доказательство неалиасинга через
// __restrict__.
//
// Стратегия: СОСУЩЕСТВОВАНИЕ AoS+SoA.
// ParticleSystem живёт всю симуляцию рядом с WorkBuffers (в
// time_integration). Перед расчётом потребителей par вызывается gather()
// (AoS->SoA), горячие циклы работают на плоских массивах с __restrict__,
// после — scatter() (SoA->AoS) для изменённых полей. Так миграция выполняется
// по одному модулю с проверкой корректности (1e-12) после каждого.
//
// Размерность массивов x[dim]/vx[dim]: во ВСЕХ потребителях пар обращения к
// координатам/скоростям идут только через индекс d<dim (циклы for(d<dim)) либо
// через константу 0 (ветка axis_sim). Явных vx[1]/vx[2] вне зависимости от dim
// нет, поэтому C-массив размером dim корректен и компилируется при любом dim,
// не расходуя память на неактивные компоненты при dim<3.
//
// Поля зеркалируют struct Particle (param.hpp): кинематика [dim], термодинамика,
// тензор напряжений (s*), тензор скоростей деформаций (e*), тензор поворота (r*),
// классификация.
#pragma once
#include <vector>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "param.hpp"

struct ParticleSystem
{
    // --- кинематика (по компонентам; все доступы только [0,dim)) ---
    std::vector<double> x[dim];
    std::vector<double> vx[dim];

    // --- термодинамика / скаляры ---
    std::vector<double> mass, rho, rho_0, p, u, c, s, e, hsml, dp, e_pl, D;

    // --- тензор напряжений (девиатор) ---
    std::vector<double> sxx, syy, szz, sxy, syz, sxz;

    // --- тензор скоростей деформаций ---
    std::vector<double> exx, eyy, ezz, exy, exz, eyz;

    // --- тензор скоростей поворота ---
    std::vector<double> rxy, rxz, ryz;

    // --- классификация ---
    std::vector<int> itype;
    std::vector<char>   axis;     // bool-семантика, плотное хранение

    // --- номер тела (object_id, fiber/PLAN.md §3.1) ---
    // Для k-файла — PID (parse_k_file); для встроенных генераторов (shock_tube
    // и др.) — совпадает с itype, чтобы критерий Solid-Solid-Contact на
    // BC2/BC6 ("разные материалы") не изменился (см. fiber/SPH.md, раздел
    // «Solid -- Solid»). Сама фильтрация контакта по body[] реализуется в
    // Ф4 — здесь только поле и его заполнение.
    std::vector<int> body;

    size_t n = 0;                 // актуальное число частиц

    // -----------------------------------------------------------------------
    // resize — резервирование всех массивов под new_n. Растёт только вверх;
    // после первого вызова на шагах с тем же n — no-op (size уже равен new_n).
    // -----------------------------------------------------------------------
    void resize(size_t new_n)
    {
        n = new_n;
        for (size_t d = 0; d < dim; ++d)
        {
            x[d].resize(n);
            vx[d].resize(n);
        }
        mass.resize(n);   rho.resize(n);   rho_0.resize(n);
        p.resize(n);      u.resize(n);     c.resize(n);     s.resize(n);     e.resize(n);
        hsml.resize(n);   dp.resize(n);    e_pl.resize(n);  D.resize(n);
        sxx.resize(n);    syy.resize(n);   szz.resize(n);
        sxy.resize(n);    syz.resize(n);   sxz.resize(n);
        exx.resize(n);    eyy.resize(n);   ezz.resize(n);
        exy.resize(n);    exz.resize(n);   eyz.resize(n);
        rxy.resize(n);    rxz.resize(n);   ryz.resize(n);
        itype.resize(n);
        axis.resize(n);
        body.resize(n);
    }

    size_t size() const { return n; }

    // -----------------------------------------------------------------------
    // Мост gather/scatter AoS<->SoA УПРАЗДНЁН (Этап 4.3.3): sys — первичное
    // хранилище, vector<Particle> не используется в production. Particle
    // остаётся только как тип одной частицы (test_particle_system, методы ниже).
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Методы-эквиваленты Particle::* для работы sys как ПЕРВИЧНОГО хранилища
    // В SoA «одна частица по индексу» = набор
    // значений в плоских массивах; эти методы воспроизводят семантику
    // Particle::obnul/inzializ/vuvod/inzializ_2d_symm по индексу.
    // -----------------------------------------------------------------------

    // Обнулить все поля частицы i (эквивалент Particle::obnul).
    void obnul_index(size_t i)
    {
        for (size_t d = 0; d < dim; ++d) { x[d][i] = 0.0; vx[d][i] = 0.0; }
        itype[i] = 0; mass[i] = 0.0; rho[i] = 0.0; rho_0[i] = 0.0;
        p[i] = 0.0; u[i] = 0.0; c[i] = 0.0; s[i] = 0.0; e[i] = 0.0;
        hsml[i] = 0.0; dp[i] = 0.0; e_pl[i] = 0.0; D[i] = 0.0;
        sxx[i] = 0.0; syy[i] = 0.0; szz[i] = 0.0;
        sxy[i] = 0.0; syz[i] = 0.0; sxz[i] = 0.0;
        exx[i] = 0.0; eyy[i] = 0.0; ezz[i] = 0.0;
        exy[i] = 0.0; exz[i] = 0.0; eyz[i] = 0.0;
        rxy[i] = 0.0; rxz[i] = 0.0; ryz[i] = 0.0;
        axis[i] = 0;
        body[i] = 0;
    }

    // Инициализация частицы i (эквивалент Particle::inzializ). В исходном AoS
    // inzializ делал push_back для x/vx; в SoA массивы уже размерены, поэтому
    // синоним obnul_index (обнулить все поля по индексу).
    void inzializ_index(size_t i) { obnul_index(i); }

    // Вывод частицы i в поток (эквивалент Particle::vuvod).
    void vuvod_index(size_t i, std::ofstream& filename)
    {
        filename << " " << itype[i];
        filename << " " << mass[i];
        for (size_t d = 0; d < dim; ++d) filename << " " << x[d][i];
        for (size_t d = 0; d < dim; ++d) filename << " " << vx[d][i];
        double sigma_i = std::sqrt((sxx[i] - syy[i]) * (sxx[i] - syy[i]) +
                                   (sxx[i] - szz[i]) * (sxx[i] - szz[i]) +
                                   (syy[i] - szz[i]) * (syy[i] - szz[i]) +
                                   6 * (sxy[i] * sxy[i] + syz[i] * syz[i] + sxz[i] * sxz[i])) / std::sqrt(2.0);
        filename << " " << rho[i] << " " << p[i] << " " << u[i] << " " << c[i] << " " << sigma_i
                 << " " << e_pl[i] << " " << D[i] << " " << hsml[i]
                 << " " << sxx[i] << " " << syy[i] << " " << szz[i]
                 << " " << sxy[i] << " " << syz[i] << " " << sxz[i] << "\n";
    }

    // Осесимметричное зеркалирование: частица dst = симметричный образ src
    // (эквивалент Particle::inzializ_2d_symm(buf)). Используется в virt_part.
    void inzializ_2d_symm_index(size_t dst, size_t src)
    {
        for (size_t k = 0; k < dim; ++k)
        {
            if (k == 0) { x[k][dst] = -x[k][src]; vx[k][dst] = -vx[k][src]; }
            else        { x[k][dst] =  x[k][src]; vx[k][dst] =  vx[k][src]; }
        }
        sxy[dst] = -sxy[src]; exy[dst] = exy[src]; rxy[dst] = rxy[src];
        itype[dst] = itype[src]; mass[dst] = mass[src];
        rho[dst] = rho[src]; rho_0[dst] = rho_0[src];
        p[dst] = p[src]; u[dst] = u[src]; c[dst] = c[src];
        s[dst] = s[src]; e[dst] = e[src];
        sxx[dst] = sxx[src]; syy[dst] = syy[src]; szz[dst] = szz[src];
        syz[dst] = syz[src]; sxz[dst] = sxz[src];
        exx[dst] = exx[src]; eyy[dst] = eyy[src]; ezz[dst] = ezz[src];
        eyz[dst] = eyz[src]; exz[dst] = exz[src];
        rxz[dst] = rxz[src]; ryz[dst] = ryz[src];
        D[dst] = D[src]; dp[dst] = dp[src]; hsml[dst] = hsml[src]; e_pl[dst] = e_pl[src];
        body[dst] = body[src];
    }
};
