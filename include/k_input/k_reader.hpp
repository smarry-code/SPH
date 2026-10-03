// k_reader.hpp
// Слой чтения LS-DYNA keyword (.k) и сырые данные KDeck (этап K2,
// readk/k-import-plan.md, разд. 3 Р7-Р9, 4.1). Без физики: только то, что
// записано в файле, в порядке появления (с учётом *INCLUDE).
//
// Ширины карточек - в k_cards.hpp. Любая ошибка - std::runtime_error с
// текстом "<файл>:<строка>: сообщение" (строка 1-based).
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace k_input
{

// Источник записи: индекс в KDeck::files и номер строки (1-based).
struct SrcLoc
{
    std::uint32_t file = 0;
    std::uint32_t line = 0;
};

struct KNode  { long long nid = 0; double x = 0, y = 0, z = 0; SrcLoc src; };
struct KSph   { long long nid = 0, pid = 0; double mass = 0; bool is_volume = false; SrcLoc src; };
struct KBeam  { long long eid = 0, pid = 0, n1 = 0, n2 = 0, n3 = 0; SrcLoc src; };
struct KSolid
{
    long long eid = 0, pid = 0;
    int nn = 0;                              // 4, 8 или 10 узлов
    std::array<long long, 10> n{};
    SrcLoc src;
};
struct KPart  { long long pid = 0; std::string title; SrcLoc src; };
struct KSpcNode { long long nid = 0;  int dofx = 0, dofy = 0, dofz = 0; SrcLoc src; };
struct KSpcSet  { long long nsid = 0; int dofx = 0, dofy = 0, dofz = 0; SrcLoc src; };
struct KNodeSet { std::vector<long long> nids; SrcLoc src; };  // src - первое появление SID

struct KDeck
{
    std::string title;                                  // *TITLE
    std::vector<std::string> files;                     // files[SrcLoc::file]
    std::vector<KNode>   nodes;
    std::unordered_map<long long, std::size_t> node_index;  // NID -> индекс в nodes (первый)
    std::vector<KSph>    sph;                           // NEND уже развёрнут
    std::vector<KBeam>   beams;
    std::vector<KSolid>  solids;
    std::vector<KPart>   parts;
    std::vector<KSpcNode> spc_nodes;
    std::vector<KSpcSet>  spc_sets;
    std::map<long long, KNodeSet> node_sets;            // SID -> NID (повтор SID дописывается)
    std::map<std::string, std::size_t> skipped;         // имя ключевого слова -> число блоков

    // "файл:строка" для сообщений.
    std::string where(const SrcLoc& s) const;
    const KNode* find_node(long long nid) const;
};

// Читает файл (и include) БЕЗ проверки ссылок.
KDeck read_k_deck_raw(const std::string& path);
// Читает файл и проверяет ссылки (Р8). Основная точка входа.
KDeck read_k_deck(const std::string& path);
// Разбор из строки; virtual_name - имя в сообщениях и база для относительных
// *INCLUDE (его каталог). Include-файлы читаются с диска.
KDeck read_k_deck_from_string(const std::string& text, const std::string& virtual_name,
                              bool validate = true);
// Проверка ссылок и дубликатов (Р8); бросает runtime_error с файлом и строкой.
void validate_k_deck(const KDeck& deck);

} // namespace k_input
