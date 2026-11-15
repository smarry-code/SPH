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
//   9 - PMMA, 10 - B4C, 11 - Al2O3,
//   12 - сталь 4340, 13 - сталь 1006, 14 - Al 7039, 15 - Al 5083-H116
//   (модель Джонсона-Кука, Materials::jc; индекс в jc = itype - n_JC).
// Индексация массивов металлов/полимеров ведётся через смещение n_Me (itype - n_Me).
#include <vector>
#include <cstddef>
#include <optional>
#include <string>
#include <cmath>
#include "eos_types.hpp"
#include "strength/johnson_cook.hpp"

// Тип УРС материала Джонсона-Кука.
enum class EosKind { Linear, Shock };

// Материал Джонсона-Кука (itype 12-15): УРС + прочность + разрушение + тепло.
// Единицы СИ (Па, К, кг/м^3, Дж/(кг К)).
struct JcMaterial
{
    std::string name;
    double rho0;
    EosKind eos;
    LinearEos lin;   // используется при eos == Linear
    ShockEos shock;  // используется при eos == Shock
    jc::JcStrength s;
    jc::JcDamage d;
    jc::JcThermal th;
    // Нижняя граница давления, Па (аналог Autodyn Hydro Pmin; отрицательная = растяжение).
    // eos_jc ограничивает p >= p_min; значение <= -1e20 -- без ограничения.
    double p_min = -1.01e20;
};

// Библиотека материалов JC (fiber/jc-implementation-plan.md §3.2).
// Источники:
//   - значения Autodyn: скриншоты Mary от 27.10.2026 (G, A, B, n, C, m, T_m,
//     T_r, c_p, K; на тех скриншотах разрушение у всех None, D1-D5 = 1.01e20);
//   - STEEL 4340 целиком, включая Failure JC: fiber/steel-4340-autodyn-constants.md
//     (скриншоты Mary 15.11.2026; приоритет над литературой);
//   - Johnson & Cook 1983 (плотность Al 7039, T_r и c_p стали 1006);
//   - Johnson & Cook 1985 (D1-D5 стали 4340, совпадают с Autodyn);
//   - Grytten et al. 2006 (D1-D5 Al 5083-H116, критерий Борвика);
//   - параметры Shock для 1006 и 7039 -- ориентир из литературы (Steinberg).
// Общее: edot0 = 1 1/с, форма скорости 1st Order, скорость по полной
// деформации, нагрев через энергию u, chi = 0.9.
inline std::vector<JcMaterial> make_jc_library()
{
    using jc::DamageRateForm;
    using jc::Heating;
    using jc::JcDamage;
    using jc::JcStrength;
    using jc::JcThermal;
    using jc::RateForm;
    using jc::RateSource;
    constexpr double kOff = 1.01e20; // плейсхолдер Autodyn "не задано"

    auto strength = [](double G, double A, double B, double n, double C, double m, double Tm, double Tr)
    {
        JcStrength s;
        s.G = G; s.A = A; s.B = B; s.n = n; s.C = C; s.m = m; s.Tm = Tm; s.Tr = Tr;
        s.edot0 = 1.0;
        s.rate_form = RateForm::FirstOrder;
        s.rate_source = RateSource::Total;
        return s;
    };
    auto thermal = [](double cp)
    {
        JcThermal t;
        t.cp = cp; t.heating = Heating::Energy; t.chi = 0.9;
        return t;
    };
    auto damage_on = [](double D1, double D2, double D3, double D4, double D5, double Tm_f, DamageRateForm f)
    {
        JcDamage d;
        d.enabled = true;
        d.D1 = D1; d.D2 = D2; d.D3 = D3; d.D4 = D4; d.D5 = D5;
        d.Tm_f = Tm_f; d.edot0_f = 1.0; d.form = f;
        return d;
    };
    auto damage_off = [kOff]()
    {
        JcDamage d;
        d.enabled = false;
        d.D1 = d.D2 = d.D3 = d.D4 = d.D5 = kOff;
        return d;
    };

    std::vector<JcMaterial> v;
    // 12: STEEL 4340, Linear.
    // Все константы, включая D1-D5, Tm_f = 1793 К и edot0_f = 1 1/с, сверены с библиотекой
    // Autodyn (fiber/steel-4340-autodyn-constants.md, 15.11.2026); D1-D5 совпадают с JC 1985.
    v.push_back({"STEEL 4340", 7830.0, EosKind::Linear, LinearEos{159e9}, ShockEos{},
                 strength(81.8e9, 792e6, 510e6, 0.26, 0.014, 1.03, 1793.0, 300.0),
                 damage_on(0.05, 3.44, -2.12, 0.002, 0.61, 1793.0, DamageRateForm::JC), thermal(477.0),
                 // TODO(verify): p_min = -4 ГПа -- середина окна -2..-6 ГПа, в котором E3 даёт
                 // откол в нужной плоскости (johnson-cook/03-pmin-energy.md). В библиотеке Autodyn
                 // для 4340 Pmin не задан (Linear EOS, Failure JC): заменить на литературную
                 // откольную прочность 4340 (с источником).
                 -4.0e9});
    // 13: STEEL 1006, Shock.
    // TODO(verify): T_r и c_p (JC 1983); Gamma0, c0, S1 (Steinberg) -- литературные, не проверены.
    v.push_back({"STEEL 1006", 7896.0, EosKind::Shock, LinearEos{}, ShockEos{2.17, 4569.0, 1.49},
                 strength(81.8e9, 350e6, 275e6, 0.36, 0.022, 1.00, 1811.0, 300.0), damage_off(), thermal(452.0)});
    // 14: AL 7039, Shock.
    // TODO(verify): плотность (JC 1983) и Gamma0, c0, S1 (Steinberg) -- не сверены.
    v.push_back({"AL 7039", 2770.0, EosKind::Shock, LinearEos{}, ShockEos{2.0, 5328.0, 1.338},
                 strength(27.6e9, 337e6, 343e6, 0.41, 0.010, 1.00, 877.0, 300.0), damage_off(), thermal(875.0)});
    // 15: AL5083H116, Linear, разрушение по Борвику (Grytten et al. 2006).
    v.push_back({"AL5083H116", 2700.0, EosKind::Linear, LinearEos{58.33e9}, ShockEos{},
                 strength(26.92e9, 167e6, 596e6, 0.551, 0.001, 0.859, 893.0, 293.0),
                 damage_on(0.178, 0.389, -2.25, 0.147, 0.0, 893.0, DamageRateForm::Borvik), thermal(910.0)});
    return v;
}

struct Materials
{
    // Смещение индекса материала JC: itype 12-15 -> jc[itype - n_JC].
    static constexpr int n_JC = 12;

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

    // Материалы Джонсона-Кука (itype 12-15), см. make_jc_library.
    std::vector<JcMaterial> jc = make_jc_library();
};

// itype принадлежит модели Джонсона-Кука (12-15).
inline bool is_jc_itype(int itype)
{
    return itype >= Materials::n_JC && itype < Materials::n_JC + 4;
}

// Материал JC по itype; вызывающая сторона обязана проверить is_jc_itype.
inline const JcMaterial& jc_material(const Materials& m, int itype)
{
    return m.jc.at(static_cast<size_t>(itype - Materials::n_JC));
}

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
    if (is_jc_itype(itype)) return jc_material(m, itype).s.G; // табличный G (неповреждённый)
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
//   - itype 12-15 (JC): K = lin.K (Linear) или rho0*c0^2 (Shock);
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
    else if (is_jc_itype(itype))
    {
        // Linear: K из таблицы; Shock: K = rho0 * c0^2 (наклон адиабаты при mu = 0).
        const JcMaterial& jm = jc_material(m, itype);
        K = (jm.eos == EosKind::Linear) ? jm.lin.K : jm.rho0 * jm.shock.c0 * jm.shock.c0;
        have_K = true;
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
