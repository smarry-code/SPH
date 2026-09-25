#pragma once
// materials.hpp
// Материальные константы (ADR-3 этапа 6.3 «Архитектурная чистота»).
//
// Ранее — неконстантные глобалы в param.hpp (G, sigy, aa, zz, gamm, const_JH_2),
// логически привязанные к уравнениям состояния (eos) и моделям прочности.
// Вынесены в самостоятельный модуль как ПОДОБЪЕКТ SimulationConfig:
//   - инициализируются один раз (значениями по умолчанию) и не изменяются runtime;
//   - размыкают циклическую зависимость eos.hpp <-> param.hpp (ADR-4 этапа 6.3);
//   - устраняют неконстантное глобальное состояние (требование ODR для 5.2-В).
//
// Формальная нумерация материалов (поле itype частицы):
//   3 - тантал, 4 - алюминий, 5 - LiF,
//   6 - полиуретан 160, 7 - полиуретан 320, 8 - полиуретан 626,
//   9 - PMMA, 10 - B4C, 11 - Al2O3.
// Индексация массивов металлов/полимеров ведётся через смещение n_Me (itype - n_Me).
#include <vector>
#include <cstddef>
#include <optional>

struct Materials
{
    // Смещение индекса материала: itype металла/полимера -> [itype - n_Me].
    static constexpr size_t n_Me = 3;

    // Модуль сдвига (Мизес).
    std::vector<double> G{69e9, 26e9, 5.5e+10, 5.27e7, 1.836e8, 4e+08, 3e+09};
    // Предел текучести для модели Мизеса.
    std::vector<double> sigy{450e6, 240e6, 2e+07, 1e+06, 3.0e+6, 4e+06, 7.7e+07};
    // Константы для УА Тэта и УРС Ми-Грюнайзена.
    std::vector<double> aa{1.4506e11, 1.97e+10, 1.6659e+10, 2.1953e8, 3.595e8, 3.0227e8, 1.8797e9};
    std::vector<double> zz{1.987, 4.2, 4.272, 1.245, 1.736, 3.5, 4.668};
    std::vector<double> gamm{2.6, 1.97};
    // Наборы констант модели керамики JH-2; 0 - B4C, 1 - Al2O3.
    // Порядок столбцов:
    //  [0]-плотность [1]-G [2]-P_hel [3]-s_hel [4]-T_n [5]-A [6]-N [7]-B [8]-M
    //  [9]-s_f_max [10]-C [11]-D1 [12]-D2 [13]-K1 [14]-K2 [15]-K3 [16]-T
    // ВНИМАНИЕ (строка 0, B4C): набор констант СМЕШАННЫЙ.
    //   - плотность и упругие константы (rho = 2517, G = 199.4e9, K1 = 230.5e9)
    //     взяты из статьи-источника эталона BC2: Grady & Moody, SAND96-0551,
    //     «Shock Compression Profiles in Ceramics» (1996), карбид бора
    //     производства Eagle Picher Industries — именно этот материал
    //     испытывался в опыте BC2 (оригинальный тест CE6). Прежние 2510 /
    //     197e9 / 233e9 отвечают марке Dow Chemical;
    //   - прочностные константы (P_hel, s_hel, A, N, B, M) статья НЕ приводит.
    //     Они остались прежними и откалиброваны под ПРЕЖНИЕ K/G, поэтому
    //     формально не согласованы с упругой частью. Менять их без источника
    //     нельзя — это была бы подгонка.
    std::vector<std::vector<double>> const_JH_2{
        {2517, 199.4e9, 8.71e9, 15.44e9, 0.03, 0.927, 0.67, 0.70, 0.85, 0.2, 0.005, 0.001, 0.5, 230.5e9, -593e9, 2800e9, 0.26e9},
        {3555, 1.165e11, 2.626e9, 3.561e9, 0.137, 0.930, 0.600, 0.110, 0.600, 1.0, 0.0, 0.005, 1.000, 1.692e11, 0.0, 0.0, 0.2e9}};
};

// Модуль сдвига по номеру материала (поле itype частицы).
// Реестр разложен по двум массивам; диспетчер воспроизводит логику
// int_force.cpp:636-652. Возвращает 0.0 для материалов без G
// (газ/вода, itype 1-2, виртуальные частицы itype < 0) — вызывающая
// сторона обязана трактовать 0 как «G неизвестен», а не как нулевую жёсткость.
inline double shear_modulus_by_itype(const Materials& m, int itype)
{
    if (itype >= 3 && itype <= 9)
    {
        const size_t k = static_cast<size_t>(itype) - Materials::n_Me;
        return k < m.G.size() ? m.G[k] : 0.0;
    }
    if (itype == 10 || itype == 11)
    {
        const size_t k = static_cast<size_t>(itype - 10);
        return k < m.const_JH_2.size() ? m.const_JH_2[k][1] : 0.0;
    }
    return 0.0;
}

// Упругие константы (E, nu) Solid-частицы -- fiber/SPH.md, «Fiber-Solid-
// Contact», «Постановка проекта: коэффициент K_fs»:
//   E = 9*K*G/(3*K+G),  nu = (3*K-2*G)/(2*(3*K+G)),
// где G -- shear_modulus_by_itype (уже реализован выше), K -- модуль
// объёмного сжатия:
//   - itype 3-9 (УРС Ми-Грюнайзена/Тэт, p = a*((rho/rho0)^z - 1)):
//     K = a*z. Проверено прямым выводом из eos.cpp::eos_ua_teta: при
//     rho=rho0 (w=1) c^2 = a*z/rho0, значит K = rho0*c^2 = a*z (формула
//     p=a*((rho/rho0)^z-1) -- частный случай УРС Тэта; для Ми-Грюнайзена
//     (eos_mi_grunaizen) та же величина a*z играет роль "холодного" модуля
//     сжатия при u=0, rho=rho0 -- используется единая формула K=a*z для
//     обоих УРС, как того требует fiber/SPH.md);
//   - itype 10-11 (керамика JH-2): K = K1 = const_JH_2[itype-10][13]
//     (см. таблицу столбцов в комментарии у Materials::const_JH_2 выше);
//   - itype 1-2 (газ/вода) и любой прочий itype без G или K -- НЕ ОПРЕДЕЛЕНО
//     (std::nullopt). Fiber-Solid-Contact с такими частицами не формируется
//     (fiber/SPH.md, тот же раздел, последний абзац подраздела K_fs).
//
// При повреждении D>0 конкретной частицы используются НЕПОВРЕЖДЁННЫЕ
// константы: функция принимает только itype (индекс в таблице материалов),
// не читает состояние (D) частицы -- K и G берутся из таблицы материала,
// поэтому результат по построению не зависит от текущего повреждения.
struct ElasticProps
{
    double E = 0.0;
    double nu = 0.0;
};

inline std::optional<ElasticProps> elastic_by_itype(const Materials& m, int itype)
{
    double K = 0.0;
    bool have_K = false;
    if (itype >= 3 && itype <= 9)
    {
        const size_t k = static_cast<size_t>(itype) - Materials::n_Me;
        if (k < m.aa.size() && k < m.zz.size())
        {
            K = m.aa[k] * m.zz[k];
            have_K = true;
        }
    }
    else if (itype == 10 || itype == 11)
    {
        const size_t k = static_cast<size_t>(itype - 10);
        if (k < m.const_JH_2.size() && m.const_JH_2[k].size() > 13)
        {
            K = m.const_JH_2[k][13];
            have_K = true;
        }
    }
    if (!have_K || !(K > 0.0)) return std::nullopt;

    const double G = shear_modulus_by_itype(m, itype);
    if (!(G > 0.0)) return std::nullopt; // газ/вода (itype 1-2) и прочие без G

    const double denom = 3.0 * K + G;
    if (!(denom > 0.0)) return std::nullopt;
    const double E = 9.0 * K * G / denom;
    if (!(E > 0.0)) return std::nullopt;
    const double nu = (3.0 * K - 2.0 * G) / (2.0 * denom);
    return ElasticProps{E, nu};
}
