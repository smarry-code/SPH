// k_cards.hpp
// Таблица описаний карточек LS-DYNA keyword (.k) для слоя чтения k_input
// (readk/k-import-plan.md, раздел 4.1; readk/ls-dyna-k-format-spec.md, разд. 3-6).
//
// ЕДИНСТВЕННОЕ место, где живут ширины полей. Смещение поля не хранится: оно
// всегда накопленная сумма ширин предыдущих полей (standard) либо
// max(ширина, 20) (long, спец. 1.3) — см. field_offset().
//
// Пункты, помеченные [K], НЕ сверены с официальным мануалом LS-DYNA
// (этап K0 пропущен по решению владельца проекта).
#pragma once

#include <cstddef>

namespace k_input
{

enum class FieldType : char
{
    I, // целое
    R, // вещественное
    A  // строка (в данных не используется, поле пропускается)
};

struct FieldDef
{
    const char* name;
    int         width;     // ширина в standard-формате
    FieldType   type;
    bool        required;  // пустое обязательное поле - ошибка
    double      def;       // умолчание для необязательного (для I - целое значение)
};

struct CardDef
{
    const FieldDef* fields;
    int             count;
};

// Ширина поля с учётом long-формата.
constexpr int field_width(const FieldDef& f, bool is_long)
{
    return (is_long && f.width < 20) ? 20 : f.width;
}

// Смещение i-го поля в строке: накопленная сумма ширин.
constexpr int field_offset(const CardDef& c, int i, bool is_long)
{
    int off = 0;
    for (int k = 0; k < i; ++k) off += field_width(c.fields[k], is_long);
    return off;
}

constexpr int kMaxCardFields = 10;

#define K_REQ true
#define K_OPT false

// *NODE: NID 8, X/Y/Z 16, TC 8, RC 8                       [K] не сверено с мануалом LS-DYNA
inline constexpr FieldDef kNodeFields[] = {
    {"NID", 8,  FieldType::I, K_REQ, 0},
    {"X",   16, FieldType::R, K_OPT, 0.0},
    {"Y",   16, FieldType::R, K_OPT, 0.0},
    {"Z",   16, FieldType::R, K_OPT, 0.0},
    {"TC",  8,  FieldType::I, K_OPT, 0},
    {"RC",  8,  FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kNodeCard{kNodeFields, 6};

// *ELEMENT_SPH: NID 8, PID 8, MASS 16, NEND 10 (спец. 3.2; [M])
inline constexpr FieldDef kSphFields[] = {
    {"NID",  8,  FieldType::I, K_REQ, 0},
    {"PID",  8,  FieldType::I, K_REQ, 0},
    {"MASS", 16, FieldType::R, K_OPT, 0.0},
    {"NEND", 10, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kSphCard{kSphFields, 4};

// *ELEMENT_BEAM: 10 полей по 8 (спец. 6.5; [M])
inline constexpr FieldDef kBeamFields[] = {
    {"EID",   8, FieldType::I, K_REQ, 0},
    {"PID",   8, FieldType::I, K_REQ, 0},
    {"N1",    8, FieldType::I, K_REQ, 0},
    {"N2",    8, FieldType::I, K_REQ, 0},
    {"N3",    8, FieldType::I, K_OPT, 0},
    {"RT1",   8, FieldType::I, K_OPT, 0},
    {"RR1",   8, FieldType::I, K_OPT, 0},
    {"RT2",   8, FieldType::I, K_OPT, 0},
    {"RR2",   8, FieldType::I, K_OPT, 0},
    {"LOCAL", 8, FieldType::I, K_OPT, 2},
};
inline constexpr CardDef kBeamCard{kBeamFields, 10};

// *PART, карточка 1 - заголовок A80 (читается строкой целиком, не полями).
// *PART, карточка 2: PID SECID MID EOSID HGID GRAV ADPOPT TMID по 10   [K] не сверено с мануалом LS-DYNA
inline constexpr FieldDef kPartFields[] = {
    {"PID",    10, FieldType::I, K_REQ, 0},
    {"SECID",  10, FieldType::I, K_OPT, 0},
    {"MID",    10, FieldType::I, K_OPT, 0},
    {"EOSID",  10, FieldType::I, K_OPT, 0},
    {"HGID",   10, FieldType::I, K_OPT, 0},
    {"GRAV",   10, FieldType::I, K_OPT, 0},
    {"ADPOPT", 10, FieldType::I, K_OPT, 0},
    {"TMID",   10, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kPartCard{kPartFields, 8};

// *ELEMENT_SOLID: EID PID N1..N8 по 8 (спец. 6.1). N1..N4 обязательны,
// N5..N8 - необязательны (0 = нет узла).
inline constexpr FieldDef kSolidFields[] = {
    {"EID", 8, FieldType::I, K_REQ, 0},
    {"PID", 8, FieldType::I, K_REQ, 0},
    {"N1",  8, FieldType::I, K_REQ, 0},
    {"N2",  8, FieldType::I, K_REQ, 0},
    {"N3",  8, FieldType::I, K_REQ, 0},
    {"N4",  8, FieldType::I, K_REQ, 0},
    {"N5",  8, FieldType::I, K_OPT, 0},
    {"N6",  8, FieldType::I, K_OPT, 0},
    {"N7",  8, FieldType::I, K_OPT, 0},
    {"N8",  8, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kSolidCard{kSolidFields, 10};

// *ELEMENT_SOLID (ten nodes format): две карточки (спец. 6.2).
inline constexpr FieldDef kSolid10Head[] = {
    {"EID", 8, FieldType::I, K_REQ, 0},
    {"PID", 8, FieldType::I, K_REQ, 0},
};
inline constexpr CardDef kSolid10Card1{kSolid10Head, 2};
inline constexpr FieldDef kSolid10Nodes[] = {
    {"N1",  8, FieldType::I, K_REQ, 0},
    {"N2",  8, FieldType::I, K_REQ, 0},
    {"N3",  8, FieldType::I, K_REQ, 0},
    {"N4",  8, FieldType::I, K_REQ, 0},
    {"N5",  8, FieldType::I, K_REQ, 0},
    {"N6",  8, FieldType::I, K_REQ, 0},
    {"N7",  8, FieldType::I, K_REQ, 0},
    {"N8",  8, FieldType::I, K_REQ, 0},
    {"N9",  8, FieldType::I, K_REQ, 0},
    {"N10", 8, FieldType::I, K_REQ, 0},
};
inline constexpr CardDef kSolid10Card2{kSolid10Nodes, 10};

// *BOUNDARY_SPC_NODE: NID CID DOFX DOFY DOFZ DOFRX DOFRY DOFRZ по 10   [K] не сверено с мануалом LS-DYNA
inline constexpr FieldDef kSpcNodeFields[] = {
    {"NID",   10, FieldType::I, K_REQ, 0},
    {"CID",   10, FieldType::I, K_OPT, 0},
    {"DOFX",  10, FieldType::I, K_OPT, 0},
    {"DOFY",  10, FieldType::I, K_OPT, 0},
    {"DOFZ",  10, FieldType::I, K_OPT, 0},
    {"DOFRX", 10, FieldType::I, K_OPT, 0},
    {"DOFRY", 10, FieldType::I, K_OPT, 0},
    {"DOFRZ", 10, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kSpcNodeCard{kSpcNodeFields, 8};

// *BOUNDARY_SPC_SET: NSID CID DOFX..DOFRZ по 10                       [K] не сверено с мануалом LS-DYNA
inline constexpr FieldDef kSpcSetFields[] = {
    {"NSID",  10, FieldType::I, K_REQ, 0},
    {"CID",   10, FieldType::I, K_OPT, 0},
    {"DOFX",  10, FieldType::I, K_OPT, 0},
    {"DOFY",  10, FieldType::I, K_OPT, 0},
    {"DOFZ",  10, FieldType::I, K_OPT, 0},
    {"DOFRX", 10, FieldType::I, K_OPT, 0},
    {"DOFRY", 10, FieldType::I, K_OPT, 0},
    {"DOFRZ", 10, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kSpcSetCard{kSpcSetFields, 8};

// *SET_NODE_LIST (и _TITLE, _GENERATE), карточка 1: SID DA1..DA4 SOLVER по 10 (спец. 4.1)
inline constexpr FieldDef kSetHeadFields[] = {
    {"SID",    10, FieldType::I, K_REQ, 0},
    {"DA1",    10, FieldType::R, K_OPT, 0.0},
    {"DA2",    10, FieldType::R, K_OPT, 0.0},
    {"DA3",    10, FieldType::R, K_OPT, 0.0},
    {"DA4",    10, FieldType::R, K_OPT, 0.0},
    {"SOLVER", 10, FieldType::A, K_OPT, 0.0},
};
inline constexpr CardDef kSetHeadCard{kSetHeadFields, 6};

// *SET_NODE_LIST, карточки 2.. : NID1..NID8 по 10 (пустые/0 игнорируются)
inline constexpr FieldDef kSetListFields[] = {
    {"NID1", 10, FieldType::I, K_OPT, 0}, {"NID2", 10, FieldType::I, K_OPT, 0},
    {"NID3", 10, FieldType::I, K_OPT, 0}, {"NID4", 10, FieldType::I, K_OPT, 0},
    {"NID5", 10, FieldType::I, K_OPT, 0}, {"NID6", 10, FieldType::I, K_OPT, 0},
    {"NID7", 10, FieldType::I, K_OPT, 0}, {"NID8", 10, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kSetListCard{kSetListFields, 8};

// *SET_NODE_LIST_GENERATE, карточки 2..: диапазоны B1BEG..B4END по 10, включительно
inline constexpr FieldDef kSetGenFields[] = {
    {"B1BEG", 10, FieldType::I, K_OPT, 0}, {"B1END", 10, FieldType::I, K_OPT, 0},
    {"B2BEG", 10, FieldType::I, K_OPT, 0}, {"B2END", 10, FieldType::I, K_OPT, 0},
    {"B3BEG", 10, FieldType::I, K_OPT, 0}, {"B3END", 10, FieldType::I, K_OPT, 0},
    {"B4BEG", 10, FieldType::I, K_OPT, 0}, {"B4END", 10, FieldType::I, K_OPT, 0},
};
inline constexpr CardDef kSetGenCard{kSetGenFields, 8};

#undef K_REQ
#undef K_OPT

} // namespace k_input
