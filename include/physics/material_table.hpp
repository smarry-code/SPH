#pragma once
// physics/material_table.hpp
// POD-таблицы материалов -- fiber/CUDA_PLAN.md §3.2, этап Г2:
// materials.hpp::Materials хранит константы в std::vector/
// std::vector<std::vector<...>>, что запрещено в SPH_HD-коде (план §3.2,
// "Правила кода"). MaterialTable -- фиксированного размера копия РОВНО
// того, что нужно compute_solid_contact_pair (physics/solid_contact_pair.hpp,
// часть "Solid-Contact"): модуль сдвига G по itype (используется в
// pr.g_eff -- ветка Герца типа 1).
//
// ElasticTable (часть "Fiber-Solid-Contact", этот же файл по списку §3.2:
// "MaterialTable, ElasticTable, FiberMaterialTable") -- фиксированного
// размера копия того, что нужно compute_fiber_solid_contact_geometry
// (physics/fiber_solid_geom.hpp): упругие константы (E, nu) Solid-частицы.
// materials.hpp::elastic_by_itype возвращает std::optional<ElasticProps> и
// принимает Materials -- недопустимо в SPH_HD-коде так же, как и
// shear_modulus_by_itype. FiberMaterialTable (таблица по НЕСКОЛЬКИМ
// материалам НИТЕЙ) в объём этой части Г2 не входит -- fiber::FiberMaterial
// (fiber/fiber_system.hpp) уже POD (E, nu, d, rho, eps_break), вызывающая
// сторона (interactions/fiber_bond.hpp, interactions/fiber_contact.hpp,
// interactions/fiber_solid_contact.hpp) читает нужные поля КОНКРЕТНОГО
// материала пары/primitive и передаёт их как скаляры (см. обоснование в
// physics/fiber_bond_pair.hpp) -- таблица-реестр по material_id понадобится,
// когда этот выбор перейдёт на GPU (части Г6-Г8).
//
// build_material_table()/build_elastic_table() -- НЕ SPH_HD: собирают
// таблицу из Materials (std::vector/std::optional) один раз при подготовке
// параметров RHS (как MaterialTable в целом по плану §3.1: "передаётся в
// ядра по значению или через __constant__"), не вызываются из device-кода.
//
// material_shear_modulus()/elastic_table_lookup() -- построчная копия
// диспетчеров materials.hpp::shear_modulus_by_itype/elastic_by_itype (только
// индексация/диапазон -- сама формула E,nu уже посчитана в
// build_elastic_table) над POD-таблицей вместо Materials; тот же результат
// бит в бит для тех же itype.
#include <cstddef>

#include "materials.hpp"  // Materials, ElasticProps, elastic_by_itype (сборка таблиц, не SPH_HD)
#include "physics/hd.hpp" // SPH_HD

namespace physics
{

struct MaterialTable
{
    // itype 3..9 (Materials::G) -- металлы/полимеры, модуль сдвига Мизеса.
    static constexpr std::size_t kMetalCount = 7;
    // itype 10..11 (Materials::const_JH_2[.][1]) -- керамика JH-2, столбец G.
    static constexpr std::size_t kJh2Count = 2;
    // Смещение индекса itype -> [itype - kNMe] у металлов/полимеров -- та же
    // константа, что Materials::n_Me (используется физикой напрямую, без
    // объекта Materials -- значение известно на этапе компиляции).
    static constexpr std::size_t kNMe = Materials::n_Me;
    // Число столбцов в строке JH-2 -- как у Materials::const_JH_2 (см. там
    // комментарий с раскладкой столбцов).
    static constexpr std::size_t kJh2Cols = 17;

    real_t G[kMetalCount] = {};
    real_t jh2_g[kJh2Count] = {};

    // Прочность Мизеса (physics/strength.hpp), этап Г2д: предел текучести --
    // Materials::sigy, та же индексация [itype - kNMe], что у G.
    real_t sigy[kMetalCount] = {};

    // ЭОС (physics/eos.hpp), этап Г2д: aa/zz/gamm -- Materials::aa/zz/gamm
    // (УА Тэта -- itype 5..9, aa/zz; Ми-Грюнайзен -- itype 3..4, aa/zz/gamm),
    // индексация -- [itype - kNMe], буквально как у Materials. gamm заполнен
    // только у индексов 0..1 (Materials::gamm имеет 2 элемента -- используются
    // только itype 3,4), остальные -- 0 (тот же паттерн "нет данных -> 0",
    // что у G/jh2_g выше; eos_ua_teta gamm не читает, поэтому нули у
    // индексов 2..6 недостижимы).
    real_t aa[kMetalCount] = {};
    real_t zz[kMetalCount] = {};
    real_t gamm[kMetalCount] = {};

    // jh2[mat][col] -- побитная копия строки Materials::const_JH_2[mat] (все
    // 17 столбцов, не только G) -- нужна и УРС JH-2 (K1,K2,K3,T -- eos.hpp),
    // и прочности JH-2 (P_hel,s_hel,T_n,A,N,B,M,s_f_max,C,D1,D2,K1 --
    // strength.hpp).
    real_t jh2[kJh2Count][kJh2Cols] = {};
};

inline MaterialTable build_material_table(const Materials& m)
{
    MaterialTable t;
    for (std::size_t k = 0; k < MaterialTable::kMetalCount; ++k)
        t.G[k] = (k < m.G.size()) ? m.G[k] : real_t(0.0);
    for (std::size_t k = 0; k < MaterialTable::kJh2Count; ++k)
        t.jh2_g[k] = (k < m.const_JH_2.size() && m.const_JH_2[k].size() > 1) ? m.const_JH_2[k][1] : real_t(0.0);
    for (std::size_t k = 0; k < MaterialTable::kMetalCount; ++k)
        t.sigy[k] = (k < m.sigy.size()) ? m.sigy[k] : real_t(0.0);
    for (std::size_t k = 0; k < MaterialTable::kMetalCount; ++k)
    {
        t.aa[k] = (k < m.aa.size()) ? m.aa[k] : real_t(0.0);
        t.zz[k] = (k < m.zz.size()) ? m.zz[k] : real_t(0.0);
        t.gamm[k] = (k < m.gamm.size()) ? m.gamm[k] : real_t(0.0);
    }
    for (std::size_t k = 0; k < MaterialTable::kJh2Count; ++k)
    {
        for (std::size_t col = 0; col < MaterialTable::kJh2Cols; ++col)
        {
            t.jh2[k][col] = (k < m.const_JH_2.size() && col < m.const_JH_2[k].size()) ? m.const_JH_2[k][col] : real_t(0.0);
        }
    }
    return t;
}

// Побитная копия materials.hpp::shear_modulus_by_itype над MaterialTable.
SPH_HD inline real_t material_shear_modulus(const MaterialTable& t, int itype)
{
    if (itype >= 3 && itype <= 9)
    {
        const int k = itype - static_cast<int>(Materials::n_Me);
        return (k >= 0 && static_cast<std::size_t>(k) < MaterialTable::kMetalCount) ? t.G[static_cast<std::size_t>(k)]
                                                                                     : real_t(0.0);
    }
    if (itype == 10 || itype == 11)
    {
        const int k = itype - 10;
        return (k >= 0 && static_cast<std::size_t>(k) < MaterialTable::kJh2Count) ? t.jh2_g[static_cast<std::size_t>(k)]
                                                                                    : real_t(0.0);
    }
    return real_t(0.0);
}

// ----------------------------------------------------------------------
// ElasticTable -- POD-аналог MaterialTable для (E, nu) Solid-частиц, нужных
// Fiber-Solid-Contact (см. шапку файла). Та же индексация itype, что у
// MaterialTable; "не определено" (materials.hpp::elastic_by_itype вернул
// std::nullopt -- K<=0 или G<=0, напр. газ/вода itype 1-2) кодируется
// отдельным bool-флагом has[k]/jh2_has[k], т.к. std::optional запрещён в
// SPH_HD-коде.
// ----------------------------------------------------------------------
struct ElasticTable
{
    static constexpr std::size_t kMetalCount = MaterialTable::kMetalCount;
    static constexpr std::size_t kJh2Count = MaterialTable::kJh2Count;

    real_t E[kMetalCount] = {};
    real_t nu[kMetalCount] = {};
    bool has[kMetalCount] = {};

    real_t jh2_E[kJh2Count] = {};
    real_t jh2_nu[kJh2Count] = {};
    bool jh2_has[kJh2Count] = {};
};

inline ElasticTable build_elastic_table(const Materials& m)
{
    ElasticTable t;
    for (std::size_t k = 0; k < ElasticTable::kMetalCount; ++k)
    {
        const int itype = static_cast<int>(k + Materials::n_Me);
        const auto es = elastic_by_itype(m, itype);
        t.has[k] = es.has_value();
        if (es)
        {
            t.E[k] = es->E;
            t.nu[k] = es->nu;
        }
    }
    for (std::size_t k = 0; k < ElasticTable::kJh2Count; ++k)
    {
        const int itype = static_cast<int>(k + 10);
        const auto es = elastic_by_itype(m, itype);
        t.jh2_has[k] = es.has_value();
        if (es)
        {
            t.jh2_E[k] = es->E;
            t.jh2_nu[k] = es->nu;
        }
    }
    return t;
}

// Побитная копия диспетчера materials.hpp::elastic_by_itype (только
// индексация/диапазон -- сама формула E,nu уже посчитана в
// build_elastic_table) над ElasticTable вместо Materials. Возвращает false,
// если для itype нет упругих констант (itype вне диапазона ЛИБО
// elastic_by_itype вернул nullopt для него) -- E/nu не трогает (вызывающая
// сторона обязана трактовать false как "не определено", как и std::nullopt
// исходной функции).
SPH_HD inline bool elastic_table_lookup(const ElasticTable& t, int itype, real_t& E, real_t& nu)
{
    if (itype >= 3 && itype <= 9)
    {
        const int k = itype - static_cast<int>(Materials::n_Me);
        if (k >= 0 && static_cast<std::size_t>(k) < ElasticTable::kMetalCount && t.has[static_cast<std::size_t>(k)])
        {
            E = t.E[static_cast<std::size_t>(k)];
            nu = t.nu[static_cast<std::size_t>(k)];
            return true;
        }
        return false;
    }
    if (itype == 10 || itype == 11)
    {
        const int k = itype - 10;
        if (k >= 0 && static_cast<std::size_t>(k) < ElasticTable::kJh2Count && t.jh2_has[static_cast<std::size_t>(k)])
        {
            E = t.jh2_E[static_cast<std::size_t>(k)];
            nu = t.jh2_nu[static_cast<std::size_t>(k)];
            return true;
        }
        return false;
    }
    return false;
}

} // namespace physics
