#pragma once
// interactions/fiber_normal_projection_types.hpp
// Типы стадии нормального контакта законов 1-3 как ПРОЕКЦИИ СКОРОСТИ
// (fiber/LAW3_VARIANT1_VELOCITY_PROJECTION.md, fiber/TZ_LAW3_PROJECTION.md).
// Вынесены отдельно от алгоритма (interactions/fiber_normal_projection.hpp),
// чтобы их могли включать fiber/model.hpp (статистика шага), checkpoint.hpp
// (таблица тёплого старта) и kernel'ы нормального контакта без цикла include.
//
// Все структуры -- POD без указателей и std::vector: переносятся в
// device-память как есть.
#include <cstddef>
#include <cstdint>
#include "param.hpp"
#include "rhs/types.hpp"

namespace interactions
{

// Вид контакта в ключе пары тёплого старта.
inline constexpr std::int32_t kNormalKindFiberFiber = 0;
inline constexpr std::int32_t kNormalKindFiberSolid = 1;

// ----------------------------------------------------------------------
// Ключ пары для тёплого старта: (вид, k0, k1, k2, k3), сравнение
// лексикографическое. Нить-нить: (0, a0, a1, b0, b1), канонизированный --
// если a0 > b0, стороны A и B меняются местами (импульс P -- скаляр, от
// порядка сторон не зависит). Нить-Solid: (1, a, b, s, -1).
// ----------------------------------------------------------------------
struct NormalPairKey
{
    std::int32_t kind = 0;
    long long k[4] = {-1, -1, -1, -1};
};

inline bool operator<(const NormalPairKey& x, const NormalPairKey& y)
{
    if (x.kind != y.kind) return x.kind < y.kind;
    for (int i = 0; i < 4; ++i)
        if (x.k[i] != y.k[i]) return x.k[i] < y.k[i];
    return false;
}

inline bool operator==(const NormalPairKey& x, const NormalPairKey& y)
{
    if (x.kind != y.kind) return false;
    for (int i = 0; i < 4; ++i)
        if (x.k[i] != y.k[i]) return false;
    return true;
}

inline NormalPairKey make_ff_pair_key(long long a0, long long a1, long long b0, long long b1)
{
    NormalPairKey key;
    key.kind = kNormalKindFiberFiber;
    if (a0 > b0)
    {
        key.k[0] = b0; key.k[1] = b1; key.k[2] = a0; key.k[3] = a1;
    }
    else
    {
        key.k[0] = a0; key.k[1] = a1; key.k[2] = b0; key.k[3] = b1;
    }
    return key;
}

inline NormalPairKey make_fs_pair_key(long long a, long long b, long long s)
{
    NormalPairKey key;
    key.kind = kNormalKindFiberSolid;
    key.k[0] = a; key.k[1] = b; key.k[2] = s; key.k[3] = -1;
    return key;
}

// Запись таблицы тёплого старта: ключ пары и накопленный импульс прошлого шага.
struct NormalWarmEntry
{
    NormalPairKey key;
    double P = 0.0;
};

// ----------------------------------------------------------------------
// Нормальный слот нить-нить: заполняет kernel нормального контакта
// (interactions/fiber_contact.hpp) в режиме проекции для КАЖДОГО принятого
// контакта с delta > 0 (независимо от mu). Массив сжатый: только такие
// слоты, по возрастанию индекса primitive (prim). Поля a0..n -- как у FiberFrictionContact.
// P -- накопленный нормальный импульс шага (выход стадии), N = P/h.
// ----------------------------------------------------------------------
struct FiberNormalContact
{
    bool active = false; // контакт принят и delta > 0
    rhs::index_t a0 = 0;
    long long a1 = -1; // -1 -- точечный сегмент
    rhs::index_t b0 = 0;
    long long b1 = -1;
    double s = 0.0;
    double t = 0.0;
    double n[dim] = {}; // нормаль от B к A
    double delta = 0.0; // проникание на x^n (диагностика)
    double mu = 0.0;    // коэффициент трения пары (для слота трения после стадии)
    NormalPairKey key;
    double P = 0.0;
    bool warm = false; // P^(0) взят из таблицы тёплого старта
    // Место в исходном массиве примитивов. Слоты хранятся СЖАТО (только активные, по возрастанию prim):
    // seg -- номер отрезка kDetSegments, в котором стоял примитив при обходе ядра (parallel_for_segments по
    // всем примитивам). Стадия пишет вклады в буфер seg и сворачивает суммы по отрезкам seg -- порядок
    // суммирования тот же, что при обходе всех примитивов (rhs::CpuExecutor::parallel_for_seg_groups).
    std::size_t prim = 0;
    std::uint16_t seg = 0;
};

// Нормальный слот нить-Solid (interactions/fiber_solid_contact.hpp). Поля
// a..n -- как у FiberSolidFrictionContact: n -- от нити к Solid.
struct FiberSolidNormalContact
{
    bool active = false;
    rhs::index_t a = 0;
    long long b = -1; // -1 -- точечный контакт
    rhs::index_t s = 0;
    double lambda = 0.0;
    double n[dim] = {};
    double delta = 0.0;
    double mu = 0.0;
    NormalPairKey key;
    double P = 0.0;
    bool warm = false;
    std::size_t prim = 0;   // как у FiberNormalContact
    std::uint16_t seg = 0;
};

// ----------------------------------------------------------------------
// Диагностика стадии за ПОСЛЕДНИЙ шаг (пишет time_integration в
// fiber::Model::proj_stats, читает монитор сцены). enabled = false -- режим
// проекции выключен, поля не заполняются.
// ----------------------------------------------------------------------
struct NormalProjectionStats
{
    bool enabled = false;
    std::size_t n_slots_ff = 0;  // активные нормальные слоты нить-нить
    std::size_t n_slots_fs = 0;  // нить-Solid
    double max_ck = 0.0;         // max счётчик c_k по бусинам
    double max_cs = 0.0;         // max счётчик c_s по частицам Solid
    double resid = 0.0;          // max_c (pdot_c^{n+1})^+ после K итераций, м/с
    double resid_rel = 0.0;      // max_c pdot_c^+/(W_c P_c), P_c > 0, W_c = sum w^2/m_k БЕЗ c_k (= 1/m_red): остаток в долях h N_c/m_red
    double resid_load = 0.0;     // (sum_c pdot_c^+/W~_c) / (sum_c P_c) по активным слотам: недостающий импульс в долях полной нагрузки
    double warm_frac = 0.0;      // доля активных слотов с тёплым стартом
    double W_stage = 0.0;        // работа стадии за шаг: sum m (|v|^2 - |v~|^2)/2, Дж (<= 0, см. run_normal_projection)
    int iters = 0;               // K
    double W_warm = 0.0;         // работа одного тёплого старта (до итераций и до масштабирования), Дж
    double warm_scale = 1.0;     // множитель тёплого старта: 1 -- W_warm <= 0; < 1 -- P^(0) уменьшен до W = 0
    // Остаток условия Синьорини ПОСЛЕ стадии трения (по итоговым скоростям шага: трение может дать
    // нормальную составляющую в соседних контактах). Заполняет normal_projection_post_residual;
    // resid_post_valid = false -- не считался.
    bool resid_post_valid = false;
    double resid_post = 0.0;      // max_c (pdot_c^{n+1})^+ по итоговым скоростям, м/с
    double resid_load_post = 0.0; // как resid_load, по итоговым скоростям
};

} // namespace interactions
