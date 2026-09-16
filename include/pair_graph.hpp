// pair_graph.hpp
// РћР±С‰РёР№ РёРЅС„СЂР°СЃС‚СЂСѓРєС‚СѓСЂРЅС‹Р№ РјРѕРґСѓР»СЊ: CSR-РїСЂРµРґСЃС‚Р°РІР»РµРЅРёРµ РїР°СЂ РІР·Р°РёРјРѕРґРµР№СЃС‚РІРёСЏ + Р¶Р°РґРЅР°СЏ
// РѕРєСЂР°СЃРєР° РіСЂР°С„Р° РєРѕРЅС„Р»РёРєС‚РѕРІ + СѓРЅРёС„РёС†РёСЂРѕРІР°РЅРЅС‹Р№ РѕР±С…РѕРґ РїРѕ С†РІРµС‚Р°Рј.
//
// РќР°Р·РЅР°С‡РµРЅРёРµ: СѓСЃС‚СЂР°РЅРёС‚СЊ #pragma omp atomic РІ С†РёРєР»Р°С… РїРѕ РїР°СЂР°Рј.
// РџСЂРѕР±Р»РµРјР° atomic: СЂР°Р·РЅС‹Рµ РїР°СЂС‹ k, РѕР±СЂР°Р±Р°С‚С‹РІР°РµРјС‹Рµ СЂР°Р·РЅС‹РјРё РїРѕС‚РѕРєР°РјРё РѕРґРЅРѕРІСЂРµРјРµРЅРЅРѕ,
// РјРѕРіСѓС‚ СЃСЃС‹Р»Р°С‚СЊСЃСЏ РЅР° РѕРґРЅСѓ С‡Р°СЃС‚РёС†Сѓ (РїР°СЂР° k1=(3,7) Рё k2=(3,9) РїРёС€СѓС‚ РІ dvxdt[3]).
// РћРїРµСЂР°С†РёСЏ += (В«РїСЂРѕС‡РёС‚Р°С‚СЊ-СЃР»РѕР¶РёС‚СЊ-Р·Р°РїРёСЃР°С‚СЊВ») РЅРµР°С‚РѕРјР°СЂРЅР° вЂ” СЃРµРіРѕРґРЅСЏ Р·Р°С‰РёС‰РµРЅР°
// #pragma omp atomic update, С‡С‚Рѕ СЃРѕР·РґР°С‘С‚ СѓР·РєРѕРµ РјРµСЃС‚Рѕ РїР°СЂР°Р»Р»РµР»РёР·РјР° (РєРѕРЅРєСѓСЂРµРЅС†РёСЏ
// Р·Р° cache-line, Р±Р»РѕРєРёСЂРѕРІРєР° SIMD).
//
// Р РµС€РµРЅРёРµ вЂ” РїСЂР°РІРёР»СЊРЅР°СЏ СЂР°СЃРєСЂР°СЃРєР° РіСЂР°С„Р° РєРѕРЅС„Р»РёРєС‚РѕРІ: РїР°СЂС‹ РѕРґРЅРѕРіРѕ С†РІРµС‚Р° РќР• РёРјРµСЋС‚
// РѕР±С‰РёС… С‡Р°СЃС‚РёС† => РєР°Р¶РґС‹Р№ РёРЅРґРµРєСЃ РјРѕРґРёС„РёС†РёСЂСѓРµС‚СЃСЏ СЂРѕРІРЅРѕ РѕРґРЅРѕР№ РїР°СЂРѕР№ РґР°РЅРЅРѕРіРѕ С†РІРµС‚Р°
// => РєРѕРЅРєСѓСЂРµРЅС‚РЅРѕР№ Р·Р°РїРёСЃРё РЅРµС‚ => atomic РЅРµ РЅСѓР¶РµРЅ. Р¦РёРєР» РїСЂРµРІСЂР°С‰Р°РµС‚СЃСЏ РІ
// В«С†РёРєР» РїРѕ С†РІРµС‚Р°Рј, РІ РєР°Р¶РґРѕРј вЂ” С‡РёСЃС‚С‹Р№ parallel for СЃ Р±Р°СЂСЊРµСЂРѕРј РјРµР¶РґСѓ С†РІРµС‚Р°РјРёВ».
//
// Р–РёР·РЅРµРЅРЅС‹Р№ С†РёРєР» (РїРѕ РѕР±СЂР°Р·С†Сѓ WorkBuffers): СЌРєР·РµРјРїР»СЏСЂ СЃРѕР·РґР°С‘С‚СЃСЏ РћР”РРќ СЂР°Р· РЅР° РІСЃСЋ
// СЃРёРјСѓР»СЏС†РёСЋ (РІ time_integration), РјРµС‚РѕРґ build() РІС‹Р·С‹РІР°РµС‚СЃСЏ РћР”РРќ СЂР°Р· Р·Р° С€Р°Рі
// РџРћРЎР›Р• Р·Р°РІРµСЂС€РµРЅРёСЏ РїРѕРёСЃРєР° СЃРѕСЃРµРґРµР№ (single_step) Рё РїРµСЂРµРёСЃРїРѕР»СЊР·СѓРµС‚СЃСЏ Р’РЎР•РњР
// РїРѕС‚СЂРµР±РёС‚РµР»СЏРјРё РїР°СЂ (int_force, density, art_visc, av_vel, ...).
//
// РџСЂРµРґРїРѕСЃС‹Р»РєР°: 0-Р±Р°Р·РѕРІР°СЏ РёРЅРґРµРєСЃР°С†РёСЏ РїР°СЂ (Р­С‚Р°Рї 4.1.2) вЂ” pair_i[k], pair_j[k]
// СЃРѕРґРµСЂР¶Р°С‚ РёРЅРґРµРєСЃС‹ С‡Р°СЃС‚РёС† РЅР°РїСЂСЏРјСѓСЋ, Р±РµР· +-1.
#pragma once
#include <vector>
#include <cstddef>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include <fstream>

// ============================================================================
// РЈСЃР»РѕРІРЅС‹Р№ РїСѓС‚СЊ atomic/coloring: РЅР°РєРѕРїР»РµРЅРёРµ РІ РїРѕР»СЏ С‡Р°СЃС‚РёС†.
// РџРѕСЂРѕРі coloring_niac_threshold (param.hpp) СЂРµС€Р°РµС‚, РєР°РєРѕР№ РїСѓС‚СЊ Р°РєС‚РёРІРµРЅ.
//   AccAtomic вЂ” РґР»СЏ РјР°Р»С‹С… niac: #pragma omp atomic update (Р±РµР·РѕРїР°СЃРЅРѕ РїСЂРё РіРѕРЅРєР°С…);
//   AccPlain  вЂ” РґР»СЏ coloring: РїР°СЂ РѕРґРЅРѕРіРѕ С†РІРµС‚Р° РЅРµ РєРѕРЅС„Р»РёРєС‚СѓСЋС‚, atomic РЅРµ РЅСѓР¶РµРЅ.
// РўРµР»Рѕ С†РёРєР»Р° РїРѕ РїР°СЂРµ РїРёС€РµС‚СЃСЏ РћР”РРќ СЂР°Р· (С€Р°Р±Р»РѕРЅРЅР°СЏ Р»СЏРјР±РґР°) Рё РёРЅСЃС‚Р°РЅС†РёСЂСѓРµС‚СЃСЏ РґР»СЏ
// РѕР±РѕРёС… Р°РєРєСѓРјСѓР»СЏС‚РѕСЂРѕРІ => РЅРµС‚ РґСѓР±Р»РёСЂРѕРІР°РЅРёСЏ С‚РµР»Р° Рё РЅРµС‚ РЅР°РєР»Р°РґРЅС‹С… СЂР°СЃС…РѕРґРѕРІ РІ runtime.
// ============================================================================
struct AccAtomic
{
    static inline void add(double& x, double v)
    {
#pragma omp atomic update
        x += v;
    }
};
struct AccPlain
{
    static inline void add(double& x, double v) { x += v; }
};

struct PairGraph
{
    // --- CSR РїРѕ С‡Р°СЃС‚РёС†Р°Рј (СѓСЃРєРѕСЂСЏРµС‚ РѕРєСЂР°СЃРєСѓ: РґР»СЏ С‡Р°СЃС‚РёС†С‹ v вЂ” СЃРїРёСЃРѕРє РїР°СЂ) ---
    std::vector<size_t> pair_offset;  // [ntotal+1]   вЂ” СЌРєСЃРєР»СЋР·РёРІРЅР°СЏ РїСЂРµС„РёРєСЃРЅР°СЏ СЃСѓРјРјР° СЃС‚РµРїРµРЅРµР№ С‡Р°СЃС‚РёС†
    std::vector<size_t> pair_csr;     // [2*niac]     вЂ” РёРЅРґРµРєСЃС‹ РїР°СЂ, СЃРіСЂСѓРїРїРёСЂРѕРІР°РЅРЅС‹Рµ РїРѕ С‡Р°СЃС‚РёС†Р°Рј

    // --- РѕРєСЂР°СЃРєР° РіСЂР°С„Р° РєРѕРЅС„Р»РёРєС‚РѕРІ ---
    int num_colors = 0;               // РёСЃРїРѕР»СЊР·РѕРІР°РЅРѕ С†РІРµС‚РѕРІ
    std::vector<int> color;           // [niac]       вЂ” С†РІРµС‚ РєР°Р¶РґРѕР№ РїР°СЂС‹ (0 .. num_colors-1)
    std::vector<size_t> color_offset;   // [num_colors+1] вЂ” РїСЂРµС„РёРєСЃРЅР°СЏ СЃСѓРјРјР° С‡РёСЃР»Р° РїР°СЂ РїРѕ С†РІРµС‚Р°Рј
    std::vector<size_t> pairs_by_color; // [niac]       вЂ” РёРЅРґРµРєСЃС‹ РїР°СЂ, РѕС‚СЃРѕСЂС‚РёСЂРѕРІР°РЅРЅС‹Рµ РїРѕ С†РІРµС‚Р°Рј

    // --- Р’РЅСѓС‚СЂРµРЅРЅРёРµ scratch-Р±СѓС„РµСЂС‹ (persistent-С‡Р»РµРЅС‹). Р Р°РЅРµРµ csr_cursor/color_count/ccursor
    // Р±С‹Р»Рё Р›РћРљРђР›Р¬РќР«РњР РїРµСЂРµРјРµРЅРЅС‹РјРё build() Рё Р°Р»Р»РѕС†РёСЂРѕРІР°Р»РёСЃСЊ РЅР° РєР°Р¶РґРѕРј С€Р°РіРµ
    // (O(ntotal)/O(num_colors) Р°Р»Р»РѕРєР°С†РёР№ Р·Р° С€Р°Рі). РўРµРїРµСЂСЊ вЂ” С‡Р»РµРЅС‹: resize вЂ” no-op РїРѕСЃР»Рµ
    // РїРµСЂРІРѕРіРѕ С€Р°РіР°, reset С‡РµСЂРµР· std::fill (Р±РµР· СЂРµР°Р»Р»РѕРєР°С†РёРё). pair_count СѓРґР°Р»С‘РЅ РєР°Рє РјС‘СЂС‚РІС‹Р№
    // (Р·Р°РїРёСЃС‹РІР°Р»СЃСЏ, РЅРѕ РЅРµ С‡РёС‚Р°Р»СЃСЏ РїРѕС‚СЂРµР±РёС‚РµР»СЏРјРё); РµРіРѕ СЂРѕР»СЊ (СЃС‚РµРїРµРЅРё в†’ РєСѓСЂСЃРѕСЂС‹ scatter)
    // РІС‹РїРѕР»РЅСЏРµС‚ csr_cursor. ---
    std::vector<size_t> csr_cursor;   // [ntotal]   вЂ” С„Р°Р·Р° A: СЃС‚РµРїРµРЅРё С‡Р°СЃС‚РёС†, Р·Р°С‚РµРј РєСѓСЂСЃРѕСЂС‹ scatter
    std::vector<int> forbidden;       // РјР°СЂРєРµСЂРЅС‹Р№ РјР°СЃСЃРёРІ Р·Р°РїСЂРµС‰С‘РЅРЅС‹С… С†РІРµС‚РѕРІ (С„Р°Р·Р° B)
    std::vector<size_t> color_count;  // [num_colors] вЂ” СЃС‡С‘С‚С‡РёРєРё РїР°СЂ РїРѕ С†РІРµС‚Р°Рј (С„Р°Р·Р° C)
    std::vector<size_t> ccursor;      // [num_colors] вЂ” РєСѓСЂСЃРѕСЂС‹ scatter РїРѕ С†РІРµС‚Р°Рј (С„Р°Р·Р° C)

    // РЎС‚Р°С‚РёСЃС‚РёРєР° С„Р°Р·С‹ B1 (РґР»СЏ Р±РµРЅС‡РјР°СЂРєРѕРІ/РґРёР°РіРЅРѕСЃС‚РёРєРё; РЅРµ РёСЃРїРѕР»СЊР·СѓРµС‚СЃСЏ РІ РіРѕСЂСЏС‡РµРј РїСѓС‚Рё).
    long stat_conflicts = 0;          // С‡РёСЃР»Рѕ РєРѕРЅС„Р»РёРєС‚СѓСЋС‰РёС… РїР°СЂ (РїРµСЂРµРЅР°Р·РЅР°С‡РµРЅРЅС‹С… РІ С„Р°Р·Рµ c)
    long stat_phase_a_us = 0;         // РІСЂРµРјСЏ С„Р°Р·С‹ (a) СЃРїРµРєСѓР»СЏС‚РёРІРЅРѕР№ РѕРєСЂР°СЃРєРё, РјРєСЃ
    long stat_phase_b_us = 0;         // РІСЂРµРјСЏ С„Р°Р·С‹ (b) РѕР±РЅР°СЂСѓР¶РµРЅРёСЏ РєРѕРЅС„Р»РёРєС‚РѕРІ, РјРєСЃ
    long stat_phase_c_us = 0;         // РІСЂРµРјСЏ С„Р°Р·С‹ (c) СЃРµСЂРёР№РЅРѕРіРѕ СЂР°Р·СЂРµС€РµРЅРёСЏ, РјРєСЃ

    // --- Р‘СѓС„РµСЂС‹ С„Р°Р·С‹ B1 (СЃРїРµРєСѓР»СЏС‚РёРІРЅР°СЏ Р±Р»РѕС‡РЅРѕ-РїР°СЂР°Р»Р»РµР»СЊРЅР°СЏ РѕРєСЂР°СЃРєР°, Р­С‚Р°Рї 1.2) ---
    std::vector<int> conflict_flag;   // [niac] вЂ” С„Р»Р°Рі РєРѕРЅС„Р»РёРєС‚Р° РїР°СЂС‹ (С‚СЂРµР±СѓРµС‚ СЃРµСЂРёР№РЅРѕРіРѕ СЂР°Р·СЂРµС€РµРЅРёСЏ)
    std::vector<int> fb_par;          // [nthreads*fb_stride] вЂ” per-thread forbidden (СЌРїРѕС…Р° РїРѕ РёРЅРґРµРєСЃСѓ k)
    int fb_stride = 0;                // С€Р°Рі РёРЅРґРµРєСЃР°С†РёРё fb_par (>= РІРµСЂС…РЅСЏСЏ РѕС†РµРЅРєР° С‡РёСЃР»Р° С†РІРµС‚РѕРІ)

    // --- largest-degree-first РїРѕСЂСЏРґРѕРє РѕРєСЂР°СЃРєРё ---
    // pair_order[k_pos] вЂ” РёРЅРґРµРєСЃ РїР°СЂС‹; РїР°СЂС‹ СѓРїРѕСЂСЏРґРѕС‡РµРЅС‹ РїРѕ СѓР±С‹РІР°РЅРёСЋ СЃС‚РµРїРµРЅРё
    // (deg(i)+deg(j)). Р–Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР° РІ СЌС‚РѕРј РїРѕСЂСЏРґРєРµ РґР°С‘С‚ fewer colors
    // (WelshвЂ“Powell): В«С‚СЂСѓРґРЅС‹РµВ» (РІС‹СЃРѕРєР°СЏ СЃС‚РµРїРµРЅСЊ) РїР°СЂС‹ РєСЂР°СЃСЏС‚СЃСЏ РїРµСЂРІС‹РјРё, Р·Р°Р±РёСЂР°СЏ
    // РЅРёР·РєРёРµ С†РІРµС‚Р°, В«Р»С‘РіРєРёРµВ» РґРѕСѓРєСЂР°С€РёРІР°СЋС‚СЃСЏ first-fit. РњРµРЅСЊС€Рµ С†РІРµС‚РѕРІ => РјРµРЅСЊС€Рµ
    // Р±Р°СЂСЊРµСЂРѕРІ РјРµР¶РґСѓ С†РІРµС‚Р°РјРё РІ РїРѕС‚СЂРµР±РёС‚РµР»СЏС… (con_density/int_force/art_visc/av_vel).
    std::vector<size_t> pair_order;   // [niac] вЂ” РїР°СЂС‹ РІ РїРѕСЂСЏРґРєРµ СѓР±С‹РІР°РЅРёСЏ СЃС‚РµРїРµРЅРё
    std::vector<size_t> pdeg_cnt;     // РіРёСЃС‚РѕРіСЂР°РјРјР° СЃС‚РµРїРµРЅРµР№ РїР°СЂ (bucket sort)
    std::vector<size_t> pdeg_off;     // СЃРјРµС‰РµРЅРёСЏ bucket sort (РїРѕ СѓР±С‹РІР°РЅРёСЋ СЃС‚РµРїРµРЅРё)
    std::vector<size_t> pdeg_cur;     // РєСѓСЂСЃРѕСЂС‹ scatter bucket sort

    // --- Р­С‚Р°Рї 2.3 (C2): СЃР»РёСЏРЅРёРµ СЃРѕРІРјРµСЃС‚РёРјС‹С… С†РІРµС‚РѕРІ (color coalescing) ---
    // РџРѕСЃР»Рµ РѕСЃРЅРѕРІРЅРѕР№ РѕРєСЂР°СЃРєРё С†РІРµС‚Р° СЃ РЅРµРїРµСЂРµСЃРµРєР°СЋС‰РёРјРёСЃСЏ РјРЅРѕР¶РµСЃС‚РІР°РјРё С‡Р°СЃС‚РёС† РјРѕР¶РЅРѕ
    // РѕР±СЉРµРґРёРЅРёС‚СЊ РІ В«РіСЂСѓРїРїСѓВ» (СЃСѓРїРµСЂ-С†РІРµС‚): РїР°СЂС‹ РіСЂСѓРїРїС‹ РїРѕ-РїСЂРµР¶РЅРµРјСѓ РЅРµ РєРѕРЅС„Р»РёРєС‚СѓСЋС‚
    // => РѕР±СЂР°Р±Р°С‚С‹РІР°СЋС‚СЃСЏ РѕРґРЅРёРј parallel for Р±РµР· Р±Р°СЂСЊРµСЂР° РјРµР¶РґСѓ РёСЃС…РѕРґРЅС‹РјРё С†РІРµС‚Р°РјРё.
    // РЎР¶РёРјР°РµС‚ num_colors => РјРµРЅСЊС€Рµ Р±Р°СЂСЊРµСЂРѕРІ РІ РїРѕС‚СЂРµР±РёС‚РµР»СЏС…. РРЅРІР°СЂРёР°РЅС‚ СЃРѕС…СЂР°РЅС‘РЅ.
    std::vector<char> color_conflict; // [num_colors*num_colors] вЂ” РјР°С‚СЂРёС†Р° РєРѕРЅС„Р»РёРєС‚РѕРІ С†РІРµС‚РѕРІ
    std::vector<int>  color_group;    // [num_colors] вЂ” РіСЂСѓРїРїР° (РЅРѕРІС‹Р№ С†РІРµС‚) РєР°Р¶РґРѕРіРѕ РёСЃС…РѕРґРЅРѕРіРѕ С†РІРµС‚Р°
    std::vector<int>  color_scratch;  // Р±СѓС„РµСЂ С†РІРµС‚РѕРІ РїР°СЂ РїРѕРґ С‚РµРєСѓС‰РµР№ С‡Р°СЃС‚РёС†РµР№

    // --- Р­С‚Р°Рї 2.2 (C3): РґРёР°РіРЅРѕСЃС‚РёРєР° РЅРµСЂР°РІРЅРѕРјРµСЂРЅРѕСЃС‚Рё СЂР°Р·РјРµСЂРѕРІ С†РІРµС‚РѕРІ ---
    // РџСЂРё СЃРёР»СЊРЅРѕРј РґРёСЃР±Р°Р»Р°РЅСЃРµ (max >> min) Р±Р°СЂСЊРµСЂС‹ РјРµР¶РґСѓ С†РІРµС‚Р°РјРё РїСЂРѕСЃС‚Р°РёРІР°СЋС‚;
    // РїСЂРµРґРІС‹С‡РёСЃР»РµРЅРЅС‹Рµ РіСЂР°РЅРёС†С‹ РїРѕС‚РѕРєРѕРІ (color_chunk_range) РґР°СЋС‚ РёРґРµР°Р»СЊРЅРѕ СЂР°РІРЅС‹Р№
    // static-split РїРѕ С‡РёСЃР»Сѓ РїР°СЂ (РґР»СЏ СЂР°РІРЅРѕРІРµСЃРЅС‹С… РїР°СЂ СЌРєРІРёРІР°Р»РµРЅС‚ schedule(static)).
    size_t min_color_size = 0;        // РјРёРЅ. С‡РёСЃР»Рѕ РїР°СЂ РІ С†РІРµС‚Рµ (РґРёР°РіРЅРѕСЃС‚РёРєР° РґРёСЃР±Р°Р»Р°РЅСЃР°)
    size_t max_color_size = 0;        // РјР°РєСЃ. С‡РёСЃР»Рѕ РїР°СЂ РІ С†РІРµС‚Рµ
    // Р”РёР°РіРЅРѕСЃС‚РёРєР° В«РґРѕ/РїРѕСЃР»Рµ coalescingВ»: РїРѕР·РІРѕР»СЏРµС‚ РїРѕРЅСЏС‚СЊ, СЃРѕР·РґР°С‘С‚ Р»Рё РґРёСЃР±Р°Р»Р°РЅСЃ
    // РёСЃС…РѕРґРЅР°СЏ РѕРєСЂР°СЃРєР° (largest-degree-first) РёР»Рё СЃР»РёСЏРЅРёРµ С†РІРµС‚РѕРІ.
    int stat_pre_coal_colors = 0;     // С‡РёСЃР»Рѕ С†РІРµС‚РѕРІ РґРѕ coalescing
    size_t stat_pre_coal_min = 0;     // РјРёРЅ. СЂР°Р·РјРµСЂ С†РІРµС‚Р° РґРѕ coalescing
    size_t stat_pre_coal_max = 0;     // РјР°РєСЃ. СЂР°Р·РјРµСЂ С†РІРµС‚Р° РґРѕ coalescing

    // -----------------------------------------------------------------------
    // ensure вЂ” СЂРµР·РµСЂРІРёСЂРѕРІР°РЅРёРµ РѕСЃРЅРѕРІРЅС‹С… Р±СѓС„РµСЂРѕРІ РїРѕРґ СѓРєР°Р·Р°РЅРЅС‹Рµ СЂР°Р·РјРµСЂС‹.
    // РђР»Р»РѕРєР°С†РёСЏ С‚РѕР»СЊРєРѕ РїСЂРё СЂРѕСЃС‚Рµ; resize вЂ” no-op, РµСЃР»Рё РІРјРµСЃС‚РёРјРѕСЃС‚СЊ РґРѕСЃС‚Р°С‚РѕС‡РЅР°.
    // Р’С‹Р·С‹РІР°РµС‚СЃСЏ РёР· build(); РјРѕР¶РµС‚ РІС‹Р·С‹РІР°С‚СЊСЃСЏ РёР·РІРЅРµ РґР»СЏ В«РїСЂРѕРіСЂРµРІР°В» РїРµСЂРµРґ С€Р°РіР°РјРё.
    // -----------------------------------------------------------------------
    void ensure(size_t ntotal, size_t niac)
    {
        pair_offset.resize(ntotal + 1);
        pair_csr.resize(2 * niac);
        color.resize(niac);
        pairs_by_color.resize(niac);
        csr_cursor.resize(ntotal);
        conflict_flag.resize(niac);
        pair_order.resize(niac);          // Р­С‚Р°Рї 2.1: РїРѕСЂСЏРґРѕРє largest-degree-first
        // forbidden, color_count, ccursor, fb_par СЂР°СЃС‚СѓС‚ РїРѕ С‚СЂРµР±РѕРІР°РЅРёСЋ РІ build().
    }

    // -----------------------------------------------------------------------
    // build вЂ” РїРѕСЃС‚СЂРѕРµРЅРёРµ CSR + Р¶Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР°. Р’С‹Р·С‹РІР°РµС‚СЃСЏ РћР”РРќ СЂР°Р· Р·Р° С€Р°Рі
    // РџРћРЎР›Р• Р·Р°РІРµСЂС€РµРЅРёСЏ РїРѕРёСЃРєР° СЃРѕСЃРµРґРµР№.
    //   pair_i, pair_j вЂ” СЃРїРёСЃРєРё РёРЅРґРµРєСЃРѕРІ С‡Р°СЃС‚РёС† РґР»СЏ РєР°Р¶РґРѕР№ РїР°СЂС‹ (0-Р±Р°Р·РѕРІС‹Рµ)
    //   niac           вЂ” С‡РёСЃР»Рѕ РїР°СЂ РІР·Р°РёРјРѕРґРµР№СЃС‚РІРёСЏ
    //   ntotal         вЂ” С‡РёСЃР»Рѕ С‡Р°СЃС‚РёС† (РІРєР»СЋС‡Р°СЏ РІРёСЂС‚СѓР°Р»СЊРЅС‹Рµ, РµСЃР»Рё РµСЃС‚СЊ)
    // РЎР»РѕР¶РЅРѕСЃС‚СЊ: O(niac + ntotal) РґР»СЏ CSR, O(niac * СЃСЂРµРґРЅСЏСЏ_СЃС‚РµРїРµРЅСЊ) РґР»СЏ РѕРєСЂР°СЃРєРё.
    // -----------------------------------------------------------------------
    void build(const std::vector<size_t>& pair_i,
               const std::vector<size_t>& pair_j,
               size_t niac, size_t ntotal)
    {
        // ---- Р’С‹СЂРѕР¶РґРµРЅРЅС‹Р№ СЃР»СѓС‡Р°Р№: РЅРµС‚ РїР°СЂ ----
        if (niac == 0)
        {
            pair_offset.assign(ntotal + 1, 0);
            pair_csr.clear();
            color.clear();
            color_offset.assign(1, 0);
            pairs_by_color.clear();
            num_colors = 0;
            return;
        }

        // --- Р РµР·РµСЂРІРёСЂРѕРІР°РЅРёРµ Р±СѓС„РµСЂРѕРІ (Р°Р»Р»РѕРєР°С†РёСЏ С‚РѕР»СЊРєРѕ РЅР° РїРµСЂРІРѕРј С€Р°РіРµ; resize вЂ” no-op РґР°Р»РµРµ) ---
        ensure(ntotal, niac);

#ifdef SPH_PROFILE
        // --- Р’РЅСѓС‚СЂРёС„Р°Р·РѕРІС‹Р№ С‚Р°Р№РјРёРЅРі build (РґРёР°РіРЅРѕСЃС‚РёРєР° СѓР·РєРѕРіРѕ РјРµСЃС‚Р°; РІРЅРµ РіРѕСЂСЏС‡РµРіРѕ РїСѓС‚Рё РїРѕ СѓРјРѕР»С‡Р°РЅРёСЋ) ---
        static double g_t_csr = 0, g_t_color = 0, g_t_sort = 0;
        static size_t g_build_n = 0;
        double _bt0 = omp_get_wtime();
#endif
        // ==== Р¤Р°Р·Р° A: CSR РїРѕ С‡Р°СЃС‚РёС†Р°Рј (counting sort) ====
        // csr_cursor СЃРЅР°С‡Р°Р»Р° С…СЂР°РЅРёС‚ СЃС‚РµРїРµРЅРё С‡Р°СЃС‚РёС†, Р·Р°С‚РµРј РїРµСЂРµР·Р°РїРёСЃС‹РІР°РµС‚СЃСЏ РєСѓСЂСЃРѕСЂР°РјРё scatter
        // (СЂР°РЅСЊС€Рµ РґР»СЏ СЌС‚РѕРіРѕ Р±С‹Р»Рё РѕС‚РґРµР»СЊРЅС‹Рµ pair_count + Р»РѕРєР°Р»СЊРЅС‹Р№ cursor вЂ” РґРІРµ Р°Р»Р»РѕРєР°С†РёРё/С€Р°Рі).
        std::fill_n(csr_cursor.begin(), ntotal, 0);
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k], j = pair_j[k];
            assert(i < ntotal && j < ntotal && "РРЅРґРµРєСЃ РїР°СЂС‹ РІРЅРµ РґРёР°РїР°Р·РѕРЅР° [0, ntotal)");
            csr_cursor[i]++;
            csr_cursor[j]++;
        }
        // Р­РєСЃРєР»СЋР·РёРІРЅР°СЏ РїСЂРµС„РёРєСЃРЅР°СЏ СЃСѓРјРјР° -> pair_offset.
        pair_offset[0] = 0;
        for (size_t v = 0; v < ntotal; ++v)
            pair_offset[v + 1] = pair_offset[v] + csr_cursor[v];
        // Scatter: РєСѓСЂСЃРѕСЂ Р·Р°РїРёСЃРё С‡Р°СЃС‚РёС†С‹ v РёРЅРёС†РёР°Р»РёР·РёСЂСѓРµС‚СЃСЏ РµС‘ СЃРјРµС‰РµРЅРёРµРј pair_offset[v].
        for (size_t v = 0; v < ntotal; ++v)
            csr_cursor[v] = pair_offset[v];
        for (size_t k = 0; k < niac; ++k)
        {
            size_t i = pair_i[k], j = pair_j[k];
            pair_csr[csr_cursor[i]++] = k;
            pair_csr[csr_cursor[j]++] = k;
        }

        // РџРѕСЂСЏРґРѕРє largest-degree-first (Р­С‚Р°Рї 2.1). РЎС‚СЂРѕРёС‚СЃСЏ РїРѕСЃР»Рµ CSR (РЅСѓР¶РµРЅ
        // pair_offset РґР»СЏ СЃС‚РµРїРµРЅРµР№). РСЃРїРѕР»СЊР·СѓРµС‚СЃСЏ color_serial; color_parallel
        // СЃРѕС…СЂР°РЅСЏРµС‚ CSR-Р»РѕРєР°Р»СЊРЅС‹Р№ РїРѕСЂСЏРґРѕРє (СЃРј. РєРѕРјРјРµРЅС‚Р°СЂРёР№ РІ color_parallel), РЅРѕ
        // pair_order РІСЃРµРіРґР° РІР°Р»РёРґРµРЅ вЂ” РґР»СЏ С‚РµСЃС‚РѕРІ, РґРёР°РіРЅРѕСЃС‚РёРєРё Рё Р±СѓРґСѓС‰РµРіРѕ warm-start.
        build_pair_order(pair_i, pair_j, niac);

#ifdef SPH_PROFILE
        g_t_csr += omp_get_wtime() - _bt0;
        double _bt1 = omp_get_wtime();
#endif
        // ==== Р¤Р°Р·Р° B: Р¶Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР° РіСЂР°С„Р° РєРѕРЅС„Р»РёРєС‚РѕРІ ====
        // Р“СЂР°С„: РІРµСЂС€РёРЅС‹ вЂ” РїР°СЂС‹; СЂРµР±СЂРѕ (k1,k2) РµСЃР»Рё РїР°СЂС‹ СЂР°Р·РґРµР»СЏСЋС‚ С‡Р°СЃС‚РёС†Сѓ.
        // Р–Р°РґРЅС‹Р№ Р°Р»РіРѕСЂРёС‚Рј: РґР»СЏ РїР°СЂС‹ k Р±РµСЂС‘Рј РЅР°РёРјРµРЅСЊС€РёР№ С†РІРµС‚, РЅРµ Р·Р°РЅСЏС‚С‹Р№
        // СѓР¶Рµ СЂР°СЃРєСЂР°С€РµРЅРЅС‹РјРё СЃРѕСЃРµРґСЏРјРё (РїР°СЂР°РјРё, СЂР°Р·РґРµР»СЏСЋС‰РёРјРё i РёР»Рё j).
        //
        // Р”РёСЃРїРµС‚С‡РµСЂРёР·Р°С†РёСЏ РїРѕ С‡РёСЃР»Сѓ РїРѕС‚РѕРєРѕРІ:
        //   - 1 РїРѕС‚РѕРє  -> color_serial (Р»РѕРєР°Р»СЊРЅРѕСЃС‚СЊ-РѕСЂРёРµРЅС‚РёСЂРѕРІР°РЅРЅС‹Р№ РѕР±С…РѕРґ РїРѕ pair_csr);
        //   - >1 РїРѕС‚РѕРєР° -> color_parallel (СЃРїРµРєСѓР»СЏС‚РёРІРЅР°СЏ Р±Р»РѕС‡РЅРѕ-РїР°СЂР°Р»Р»РµР»СЊРЅР°СЏ РѕРєСЂР°СЃРєР°
        //     GebremedhinвЂ“Manne: РїР°СЂР°Р»Р»РµР»СЊРЅС‹Р№ Р¶Р°РґРЅС‹Р№ РїСЂРѕС…РѕРґ + СЃРµСЂРёР№РЅРѕРµ СЂР°Р·СЂРµС€РµРЅРёРµ
        //     РєРѕРЅС„Р»РёРєС‚РѕРІ). Р”Р»СЏ SPH-РіСЂР°С„Р° (СЃС‚РµРїРµРЅСЊ ~4вЂ“8) РґР°С‘С‚ ~PГ— СѓСЃРєРѕСЂРµРЅРёРµ С„Р°Р·С‹ B.
        if (omp_get_max_threads() > 1)
            color_parallel(pair_i, pair_j, niac, ntotal);
        else
            color_serial(pair_i, pair_j, niac, ntotal);

#ifdef SPH_PROFILE
        g_t_color += omp_get_wtime() - _bt1;
        double _bt2 = omp_get_wtime();
#endif
        // (num_colors/ color[] Р·Р°РїРѕР»РЅРµРЅС‹ РјРµС‚РѕРґРѕРј color_serial/color_parallel РІС‹С€Рµ)
        // ==== Р¤Р°Р·Р° C: РіСЂСѓРїРїРёСЂРѕРІРєР° РїР°СЂ РїРѕ С†РІРµС‚Р°Рј (counting sort) ====
        // pairs_by_color[ color_offset[c] .. color_offset[c+1] ) вЂ” РІСЃРµ РїР°СЂС‹ С†РІРµС‚Р° c.
        // color_count/ccursor вЂ” persistent-С‡Р»РµРЅС‹ (СЂР°РЅСЊС€Рµ вЂ” Р»РѕРєР°Р»СЊРЅС‹Рµ Р°Р»Р»РѕРєР°С†РёРё РЅР° РєР°Р¶РґРѕРј С€Р°РіРµ).
        color_count.assign(static_cast<size_t>(num_colors), 0);
        for (size_t k = 0; k < niac; ++k)
            color_count[static_cast<size_t>(color[k])]++;
        color_offset.resize(static_cast<size_t>(num_colors) + 1);
        color_offset[0] = 0;
        for (int c = 0; c < num_colors; ++c)
        {
            const size_t color_index = static_cast<size_t>(c);
            color_offset[color_index + 1] = color_offset[color_index] + color_count[color_index];
        }
        ccursor.assign(static_cast<size_t>(num_colors), 0);
        for (size_t k = 0; k < niac; ++k)
        {
            int c = color[k];
            const size_t color_index = static_cast<size_t>(c);
            pairs_by_color[color_offset[color_index] + ccursor[color_index]++] = k;
        }

        // Р”РёР°РіРЅРѕСЃС‚РёРєР° СЂР°СЃРєСЂР°СЃРєРё Р”Рћ coalescing (РІС‹СЏРІРёС‚СЊ РёСЃС‚РѕС‡РЅРёРє РґРёСЃР±Р°Р»Р°РЅСЃР°).
        stat_pre_coal_colors = num_colors;
        stat_pre_coal_min = (num_colors > 0) ? color_offset[static_cast<size_t>(num_colors)] : 0;
        stat_pre_coal_max = 0;
        for (int c = 0; c < num_colors; ++c)
        {
            const size_t color_index = static_cast<size_t>(c);
            size_t sz = color_offset[color_index + 1] - color_offset[color_index];
            if (sz < stat_pre_coal_min) stat_pre_coal_min = sz;
            if (sz > stat_pre_coal_max) stat_pre_coal_max = sz;
        }

        // ==== СЃР»РёСЏРЅРёРµ СЃРѕРІРјРµСЃС‚РёРјС‹С… С†РІРµС‚РѕРІ (color coalescing) ====
        // Р¦РІРµС‚Р° СЃ РЅРµРїРµСЂРµСЃРµРєР°СЋС‰РёРјРёСЃСЏ РјРЅРѕР¶РµСЃС‚РІР°РјРё С‡Р°СЃС‚РёС† РѕР±СЉРµРґРёРЅСЏСЋС‚СЃСЏ РІ РіСЂСѓРїРїСѓ:
        // РїР°СЂС‹ РіСЂСѓРїРїС‹ РЅРµ РєРѕРЅС„Р»РёРєС‚СѓСЋС‚ => РѕРґРёРЅ parallel for Р±РµР· Р±Р°СЂСЊРµСЂР° РјРµР¶РґСѓ
        // РёСЃС…РѕРґРЅС‹РјРё С†РІРµС‚Р°РјРё. РЎР¶РёРјР°РµС‚ num_colors => РјРµРЅСЊС€Рµ Р±Р°СЂСЊРµСЂРѕРІ РІ РїРѕС‚СЂРµР±РёС‚РµР»СЏС….
        // РРЅРІР°СЂРёР°РЅС‚ is_valid_coloring СЃРѕС…СЂР°РЅСЏРµС‚СЃСЏ РїРѕ РїРѕСЃС‚СЂРѕРµРЅРёСЋ РіСЂСѓРїРї.
        coalesce_colors(pair_i, pair_j, niac, ntotal);

        // ==== Р­С‚Р°Рї 2.2 (C3): РґРёР°РіРЅРѕСЃС‚РёРєР° РЅРµСЂР°РІРЅРѕРјРµСЂРЅРѕСЃС‚Рё СЂР°Р·РјРµСЂРѕРІ С†РІРµС‚РѕРІ ====
        // РЎРёР»СЊРЅС‹Р№ РґРёСЃР±Р°Р»Р°РЅСЃ (max >> min) РѕР·РЅР°С‡Р°РµС‚ РїСЂРѕСЃС‚РѕР№ РїРѕС‚РѕРєРѕРІ РЅР° Р±Р°СЂСЊРµСЂРµ РјРµР¶РґСѓ
        // С†РІРµС‚Р°РјРё. min/max СЃРѕС…СЂР°РЅСЏСЋС‚СЃСЏ РґР»СЏ Р±РµРЅС‡РјР°СЂРєРѕРІ Рё СЂРµС€РµРЅРёСЏ Рѕ Р±Р°Р»Р°РЅСЃРёСЂРѕРІРєРµ.
        min_color_size = 0;
        max_color_size = 0;
        if (num_colors > 0)
        {
            min_color_size = color_offset[static_cast<size_t>(num_colors)];  // СЃС‚Р°СЂС‚ СЃ РјР°РєСЃРёРјСѓРјР° (= niac)
            for (int c = 0; c < num_colors; ++c)
            {
                const size_t color_index = static_cast<size_t>(c);
                size_t sz = color_offset[color_index + 1] - color_offset[color_index];
                if (sz < min_color_size) min_color_size = sz;
                if (sz > max_color_size) max_color_size = sz;
            }
        }
#ifdef SPH_PROFILE
        g_t_sort += omp_get_wtime() - _bt2;
        ++g_build_n;
        if (g_build_n % 1000 == 0)
        {
            std::ofstream _pf("force_timing.csv", std::ios::app);
            if (_pf)
                _pf << "buildphases," << g_build_n << "," << niac << ","
                    << (g_t_csr / static_cast<double>(g_build_n)) * 1e6 << ","
                    << (g_t_color / static_cast<double>(g_build_n)) * 1e6 << ","
                    << (g_t_sort / static_cast<double>(g_build_n)) * 1e6 << "\n";
        }
#endif
    }

    // -----------------------------------------------------------------------
    // build_pair_order вЂ” СѓРїРѕСЂСЏРґРѕС‡РёРІР°РЅРёРµ РїР°СЂ РїРѕ СѓР±С‹РІР°РЅРёСЋ СЃС‚РµРїРµРЅРё (largest-degree-first).
    // РЎС‚РµРїРµРЅСЊ РїР°СЂС‹ k = deg(pair_i[k]) + deg(pair_j[k]), РіРґРµ deg(v) вЂ” С‡РёСЃР»Рѕ РїР°СЂ,
    // РёРЅС†РёРґРµРЅС‚РЅС‹С… С‡Р°СЃС‚РёС†Рµ v (= pair_offset[v+1]-pair_offset[v]). РЎС‚Р°Р±РёР»СЊРЅС‹Р№ counting
    // (bucket) sort O(niac + max_degree): РїР°СЂС‹ СЂР°РІРЅРѕР№ СЃС‚РµРїРµРЅРё СЃРѕС…СЂР°РЅСЏСЋС‚ РёСЃС…РѕРґРЅС‹Р№
    // РїРѕСЂСЏРґРѕРє (РґРµС‚РµСЂРјРёРЅРёР·Рј). Р’С‹СЃРѕРєР°СЏ СЃС‚РµРїРµРЅСЊ РІРїРµСЂС‘Рґ => Р¶Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР° РґР°С‘С‚ fewer
    // colors (WelshвЂ“Powell): В«С‚СЂСѓРґРЅС‹РµВ» РїР°СЂС‹ Р·Р°Р±РёСЂР°СЋС‚ РЅРёР·РєРёРµ С†РІРµС‚Р°, В«Р»С‘РіРєРёРµВ»
    // РґРѕСѓРєСЂР°С€РёРІР°СЋС‚СЃСЏ first-fit. РњРµРЅСЊС€Рµ С†РІРµС‚РѕРІ => РјРµРЅСЊС€Рµ Р±Р°СЂСЊРµСЂРѕРІ РІ РїРѕС‚СЂРµР±РёС‚РµР»СЏС….
    // -----------------------------------------------------------------------
    void build_pair_order(const std::vector<size_t>& pair_i,
                          const std::vector<size_t>& pair_j,
                          size_t niac)
    {
        auto deg = [&](size_t v) { return pair_offset[v + 1] - pair_offset[v]; };
        size_t max_pdeg = 0;
        for (size_t k = 0; k < niac; ++k)
        {
            size_t dk = deg(pair_i[k]) + deg(pair_j[k]);
            if (dk > max_pdeg) max_pdeg = dk;
        }
        // Р“РёСЃС‚РѕРіСЂР°РјРјР° РїРѕ СЃС‚РµРїРµРЅСЏРј РїР°СЂ.
        pdeg_cnt.assign(max_pdeg + 1, 0);
        for (size_t k = 0; k < niac; ++k)
            ++pdeg_cnt[deg(pair_i[k]) + deg(pair_j[k])];
        // РЎРјРµС‰РµРЅРёСЏ РґР»СЏ РџРћ РЈР‘Р«Р’РђРќРР®: СЃС‚РµРїРµРЅСЊ max_pdeg -> РїРѕР·РёС†РёСЏ 0.
        pdeg_off.assign(max_pdeg + 2, 0);
        size_t pos = 0;
        for (size_t d = max_pdeg + 1; d-- > 0;)
        {
            pdeg_off[d] = pos;
            pos += pdeg_cnt[d];
        }
        // РЎС‚Р°Р±РёР»СЊРЅС‹Р№ scatter: РїР°СЂС‹ СЂР°РІРЅРѕР№ СЃС‚РµРїРµРЅРё РёРґСѓС‚ РІ РёСЃС…РѕРґРЅРѕРј РїРѕСЂСЏРґРєРµ.
        pdeg_cur.assign(max_pdeg + 2, 0);
        for (size_t d = 0; d <= max_pdeg; ++d)
            pdeg_cur[d] = pdeg_off[d];
        for (size_t k = 0; k < niac; ++k)
        {
            size_t dk = deg(pair_i[k]) + deg(pair_j[k]);
            size_t slot = pdeg_cur[dk]++;  // operator[] РІРѕР·РІСЂР°С‰Р°РµС‚ СЃСЃС‹Р»РєСѓ (lvalue)
            pair_order[slot] = k;
        }
    }

    // -----------------------------------------------------------------------
    // coalesce_colors вЂ” СЃР»РёСЏРЅРёРµ СЃРѕРІРјРµСЃС‚РёРјС‹С… С†РІРµС‚РѕРІ (color coalescing).
    // Р”РІР° С†РІРµС‚Р° СЃРѕРІРјРµСЃС‚РёРјС‹, РµСЃР»Рё РЅРё РѕРґРЅР° РїР°СЂР° РѕРґРЅРѕРіРѕ РЅРµ СЂР°Р·РґРµР»СЏРµС‚ С‡Р°СЃС‚РёС†Сѓ СЃ РїР°СЂРѕР№
    // РґСЂСѓРіРѕРіРѕ. РЎРѕРІРјРµСЃС‚РёРјС‹Рµ С†РІРµС‚Р° РѕР±СЉРµРґРёРЅСЏСЋС‚СЃСЏ РІ РіСЂСѓРїРїСѓ (СЃСѓРїРµСЂ-С†РІРµС‚): РїР°СЂС‹ РіСЂСѓРїРїС‹
    // РїРѕ-РїСЂРµР¶РЅРµРјСѓ РЅРµ РєРѕРЅС„Р»РёРєС‚СѓСЋС‚ => РѕР±СЂР°Р±Р°С‚С‹РІР°СЋС‚СЃСЏ РѕРґРЅРёРј parallel for Р±РµР· Р±Р°СЂСЊРµСЂР°
    // РјРµР¶РґСѓ РёСЃС…РѕРґРЅС‹РјРё С†РІРµС‚Р°РјРё. РЎР¶РёРјР°РµС‚ num_colors => РјРµРЅСЊС€Рµ Р±Р°СЂСЊРµСЂРѕРІ РІ РїРѕС‚СЂРµР±РёС‚РµР»СЏС….
    // РРЅРІР°СЂРёР°РЅС‚ is_valid_coloring СЃРѕС…СЂР°РЅСЏРµС‚СЃСЏ РїРѕ РїРѕСЃС‚СЂРѕРµРЅРёСЋ РіСЂСѓРїРї.
    // РЎР»РѕР¶РЅРѕСЃС‚СЊ: O(ntotal*dМ„ВІ) РЅР° РјР°С‚СЂРёС†Сѓ РєРѕРЅС„Р»РёРєС‚РѕРІ + O(CВІ) РЅР° Р¶Р°РґРЅСѓСЋ РіСЂСѓРїРїРёСЂРѕРІРєСѓ.
    // -----------------------------------------------------------------------
    void coalesce_colors(const std::vector<size_t>& /*pair_i*/,
                         const std::vector<size_t>& /*pair_j*/,
                         size_t niac, size_t ntotal)
    {
        if (num_colors <= 1) return;  // РЅРµС‡РµРіРѕ СЃР»РёРІР°С‚СЊ
        const size_t C = static_cast<size_t>(num_colors);

        // 1. РњР°С‚СЂРёС†Р° РєРѕРЅС„Р»РёРєС‚РѕРІ С†РІРµС‚РѕРІ: РґР»СЏ РєР°Р¶РґРѕР№ С‡Р°СЃС‚РёС†С‹ С†РІРµС‚Р° РµС‘ РїР°СЂ РїРѕРїР°СЂРЅРѕ
        //    РєРѕРЅС„Р»РёРєС‚СѓСЋС‚ (СЂР°Р·РґРµР»СЏСЋС‚ СЌС‚Сѓ С‡Р°СЃС‚РёС†Сѓ).
        color_conflict.assign(C * C, 0);
        for (size_t v = 0; v < ntotal; ++v)
        {
            size_t d0 = pair_offset[v], d1 = pair_offset[v + 1];
            if (d1 - d0 < 2) continue;  // < 2 РїР°СЂ РїРѕРґ С‡Р°СЃС‚РёС†РµР№ вЂ” РјРµР¶С†РІРµС‚РЅС‹С… РєРѕРЅС„Р»РёРєС‚РѕРІ РЅРµС‚
            color_scratch.clear();
            for (size_t p = d0; p < d1; ++p)
                color_scratch.push_back(color[pair_csr[p]]);
            for (size_t a = 0; a + 1 < color_scratch.size(); ++a)
                for (size_t b = a + 1; b < color_scratch.size(); ++b)
                {
                    int ca = color_scratch[a], cb = color_scratch[b];
                    if (ca == cb) continue;  // РѕРґРёРЅ С†РІРµС‚ вЂ” РЅРµ РјРµР¶С†РІРµС‚РЅС‹Р№ РєРѕРЅС„Р»РёРєС‚
                    color_conflict[static_cast<size_t>(ca) * C + static_cast<size_t>(cb)] = 1;
                    color_conflict[static_cast<size_t>(cb) * C + static_cast<size_t>(ca)] = 1;
                }
        }
        // 2. Р–Р°РґРЅРѕРµ РЅР°Р·РЅР°С‡РµРЅРёРµ РіСЂСѓРїРї: С†РІРµС‚ c РёРґС‘С‚ РІ РїРµСЂРІСѓСЋ СЃРѕРІРјРµСЃС‚РёРјСѓСЋ РіСЂСѓРїРїСѓ.
        color_group.assign(C, -1);
        std::vector<std::vector<int>> members;  // members[g] вЂ” С†РІРµС‚Р°, РІС…РѕРґСЏС‰РёРµ РІ РіСЂСѓРїРїСѓ g
        for (size_t c = 0; c < C; ++c)
        {
            int g = -1;
            for (size_t gi = 0; gi < members.size(); ++gi)
            {
                bool ok = true;
                for (int cp : members[gi])
                    if (color_conflict[c * C + static_cast<size_t>(cp)]) { ok = false; break; }
                if (ok) { g = static_cast<int>(gi); break; }
            }
            if (g < 0) { g = static_cast<int>(members.size()); members.emplace_back(); }
            members[static_cast<size_t>(g)].push_back(static_cast<int>(c));
            color_group[c] = g;
        }
        const int num_groups = static_cast<int>(members.size());
        if (num_groups >= static_cast<int>(C))
            return;  // СЃР»РёСЏРЅРёРµ РЅРµ СѓРјРµРЅСЊС€РёР»Рѕ С‡РёСЃР»Рѕ С†РІРµС‚РѕРІ вЂ” СЂР°СЃРєСЂР°СЃРєР° СѓР¶Рµ РїР»РѕС‚РЅР°СЏ

        // 3. РџРµСЂРµРЅСѓРјРµСЂР°С†РёСЏ С†РІРµС‚Р° РїР°СЂ РїРѕ РіСЂСѓРїРїР°Рј + РїРµСЂРµСЃС‚СЂРѕР№РєР° pairs_by_color.
        for (size_t k = 0; k < niac; ++k)
            color[k] = color_group[static_cast<size_t>(color[k])];
        num_colors = num_groups;
        color_count.assign(static_cast<size_t>(num_colors), 0);
        for (size_t k = 0; k < niac; ++k)
            ++color_count[static_cast<size_t>(color[k])];
        color_offset.resize(static_cast<size_t>(num_colors) + 1);
        color_offset[0] = 0;
        for (int c = 0; c < num_colors; ++c)
        {
            const size_t color_index = static_cast<size_t>(c);
            color_offset[color_index + 1] = color_offset[color_index] + color_count[color_index];
        }
        ccursor.assign(static_cast<size_t>(num_colors), 0);
        for (size_t k = 0; k < niac; ++k)
        {
            int c = color[k];
            const size_t color_index = static_cast<size_t>(c);
            pairs_by_color[color_offset[color_index] + ccursor[color_index]++] = k;
        }
    }

    // -----------------------------------------------------------------------
    // color_chunk_range вЂ” СЃР±Р°Р»Р°РЅСЃРёСЂРѕРІР°РЅРЅС‹Р№ static-split С†РІРµС‚Р° РїРѕ РїРѕС‚РѕРєР°Рј (C3).
    // Р’РѕР·РІСЂР°С‰Р°РµС‚ [begin,end) РґРёР°РїР°Р·РѕРЅ pairs_by_color РґР»СЏ РїРѕС‚РѕРєР° tid РёР· nthreads
    // РІ С†РІРµС‚Рµ c. Р”РµР»РµРЅРёРµ СЂРѕРІРЅРѕ РїРѕ С‡РёСЃР»Сѓ РїР°СЂ (РѕСЃС‚Р°С‚РѕРє СЂР°СЃРїСЂРµРґРµР»С‘РЅ РЅР° РїРµСЂРІС‹Рµ РїРѕС‚РѕРєРё).
    // Р”Р»СЏ СЂР°РІРЅРѕРІРµСЃРЅС‹С… РїР°СЂ СЌРєРІРёРІР°Р»РµРЅС‚РЅРѕ schedule(static); С…СЂР°РЅРёС‚СЃСЏ РєР°Рє helper РґР»СЏ
    // РґРёР°РіРЅРѕСЃС‚РёРєРё РґРёСЃР±Р°Р»Р°РЅСЃР° Рё Р±СѓРґСѓС‰РµРіРѕ weighted-balancing РїСЂРё РЅРµСЂР°РІРЅРѕРјРµСЂРЅС‹С… РІРµСЃР°С….
    // -----------------------------------------------------------------------
    void color_chunk_range(size_t c, size_t tid, size_t nthreads, size_t& begin, size_t& end) const
    {
        size_t b = color_offset[c];
        size_t e = color_offset[c + 1];
        size_t W = e - b;
        size_t base = W / nthreads;
        size_t rem = W % nthreads;
        // РїРµСЂРІС‹Рµ 'rem' РїРѕС‚РѕРєРѕРІ РїРѕР»СѓС‡Р°СЋС‚ РїРѕ base+1 РїР°СЂ, РѕСЃС‚Р°Р»СЊРЅС‹Рµ вЂ” РїРѕ base.
        begin = b + tid * base + (tid < rem ? tid : rem);
        end = begin + base + (tid < rem ? 1 : 0);
    }

    // -----------------------------------------------------------------------
    // color_serial вЂ” СЃРµСЂРёР№РЅР°СЏ Р¶Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР° РІ РїРѕСЂСЏРґРєРµ largest-degree-first.
    // РћР±С…РѕРґ РїРѕ pair_order (РїР°СЂС‹ РїРѕ СѓР±С‹РІР°РЅРёСЋ СЃС‚РµРїРµРЅРё); mark_neighbors СЃРєР°РЅРёСЂСѓРµС‚
    // pair_csr РїРѕРґ С‡Р°СЃС‚РёС†Р°РјРё РїР°СЂС‹ (CSR-Р»РѕРєР°Р»СЊРЅРѕСЃС‚СЊ С‡С‚РµРЅРёСЏ СЃРѕСЃРµРґРµР№). РљР°Р¶РґР°СЏ РїР°СЂР°
    // РєСЂР°СЃРёС‚СЃСЏ СЂРѕРІРЅРѕ РѕРґРёРЅ СЂР°Р· => РїСЂРѕС‰Рµ Рё Р±С‹СЃС‚СЂРµРµ РїСЂРµР¶РЅРµРіРѕ РґРІРѕР№РЅРѕРіРѕ РѕР±С…РѕРґР° CSR.
    // -----------------------------------------------------------------------
    void color_serial(const std::vector<size_t>& pair_i,
                      const std::vector<size_t>& pair_j,
                      size_t niac, size_t /*ntotal*/)
    {
        std::fill_n(color.begin(), niac, -1);
        // РњР°СЂРєРµСЂРЅС‹Р№ РјР°СЃСЃРёРІ Р·Р°РїСЂРµС‰С‘РЅРЅС‹С… С†РІРµС‚РѕРІ: forbidden[c]==k => С†РІРµС‚ c Р·Р°РЅСЏС‚ СЃРѕСЃРµРґРѕРј k.
        std::fill(forbidden.begin(), forbidden.end(), -1);
        num_colors = 0;

        auto mark_neighbors = [&](size_t v, int cur_k) {
            for (size_t p = pair_offset[v]; p < pair_offset[v + 1]; ++p)
            {
                size_t k2 = pair_csr[p];
                if (k2 != static_cast<size_t>(cur_k) && color[k2] != -1)
                    forbidden[static_cast<size_t>(color[k2])] = cur_k;  // РїРѕРјРµС‚РёС‚СЊ С†РІРµС‚ СЃРѕСЃРµРґР° РєР°Рє Р·Р°РЅСЏС‚С‹Р№
            }
        };

        for (size_t idx = 0; idx < niac; ++idx)
        {
            size_t k = pair_order[idx];
            size_t i = pair_i[k], j = pair_j[k];
            if (static_cast<int>(forbidden.size()) < num_colors + 1)
                forbidden.resize(static_cast<size_t>(num_colors) + 1, -1);

            const int marker = static_cast<int>(k);
            mark_neighbors(i, marker);
            mark_neighbors(j, marker);

            int c = 0;
            while (c < num_colors && forbidden[static_cast<size_t>(c)] == marker)
                ++c;
            if (c == num_colors)
            {
                num_colors++;
                forbidden.push_back(-1);
            }
            color[k] = c;
        }
    }

    // -----------------------------------------------------------------------
    // color_parallel вЂ” СЃРїРµРєСѓР»СЏС‚РёРІРЅР°СЏ Р±Р»РѕС‡РЅРѕ-РїР°СЂР°Р»Р»РµР»СЊРЅР°СЏ Р¶Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР°
    // (GebremedhinвЂ“Manne 2000). РўСЂРё РїРѕРґС„Р°Р·С‹:
    //   (a) РїР°СЂР°Р»Р»РµР»СЊРЅС‹Р№ Р¶Р°РґРЅС‹Р№ РїСЂРѕС…РѕРґ РїРѕ Р±Р»РѕРєР°Рј РїР°СЂ: РєР°Р¶РґС‹Р№ РїРѕС‚РѕРє РєСЂР°СЃРёС‚ СЃРІРѕР№ Р±Р»РѕРє,
    //       РѕРїСЂРѕРјРµС‚С‡РёРІРѕ С‡РёС‚Р°СЏ color[] СЃРѕСЃРµРґРµР№ (РЅРµРєРѕС‚РѕСЂС‹Рµ РјРѕРіСѓС‚ Р±С‹С‚СЊ РµС‰С‘ РЅРµ РЅР°Р·РЅР°С‡РµРЅС‹);
    //   (b) РїР°СЂР°Р»Р»РµР»СЊРЅРѕРµ РѕР±РЅР°СЂСѓР¶РµРЅРёРµ РєРѕРЅС„Р»РёРєС‚РѕРІ (РїР°СЂС‹ РѕРґРЅРѕРіРѕ С†РІРµС‚Р°, СЂР°Р·РґРµР»СЏСЋС‰РёРµ С‡Р°СЃС‚РёС†Сѓ,
    //       Р»РёР±Рѕ РїР°СЂС‹ СЃ РїСЂРµРІС‹С€РµРЅРЅРѕР№ РѕС†РµРЅРєРѕР№ color[k] == -2);
    //   (c) СЃРµСЂРёР№РЅРѕРµ СЂР°Р·СЂРµС€РµРЅРёРµ РєРѕРЅС„Р»РёРєС‚СѓСЋС‰РёС… РїР°СЂ РЅР° СЃРІРѕР±РѕРґРЅС‹Р№ С†РІРµС‚ (СЃ СѓС‡С‘С‚РѕРј Р’РЎР•РҐ
    //       СЃРѕСЃРµРґРµР№ => РіР°СЂР°РЅС‚РёСЏ РѕС‚СЃСѓС‚СЃС‚РІРёСЏ РЅРѕРІС‹С… РєРѕРЅС„Р»РёРєС‚РѕРІ).
    // Р”Р»СЏ SPH-РіСЂР°С„Р° (СЃС‚РµРїРµРЅСЊ ~4вЂ“8) РґРѕР»СЏ РєРѕРЅС„Р»РёРєС‚РѕРІ РјР°Р»Р° => (a) РґРѕРјРёРЅРёСЂСѓРµС‚ Рё РёРґС‘С‚ РІ
    // P РїРѕС‚РѕРєРѕРІ; (c) Р·Р°С‚СЂР°РіРёРІР°РµС‚ Р»РёС€СЊ РіРѕСЂСЃС‚РєСѓ РїР°СЂ. Р§С‚РµРЅРёРµ color[] СЃРѕСЃРµРґРµР№ РЅР° x86
    // Р°С‚РѕРјР°СЂРЅРѕ (РІС‹СЂРѕРІРЅРµРЅРЅС‹Р№ int) => СЃРїРµРєСѓР»СЏС†РёСЏ Р±РµР·РѕРїР°СЃРЅР° (С…СѓРґС€РёР№ СЃР»СѓС‡Р°Р№ вЂ” Р»РёС€РЅРёР№ С†РІРµС‚,
    // СѓСЃС‚СЂР°РЅСЏРµРјС‹Р№ РІ (c)). Р—Р°РїРёСЃСЊ color[k] Р±РµСЃРєРѕРЅС„Р»РёРєС‚РЅР° (СЃС‚Р°С‚РёС‡РµСЃРєРѕРµ СЂР°СЃРїРёСЃР°РЅРёРµ:
    // СЂР°Р·РЅС‹Рµ РїРѕС‚РѕРєРё РїРёС€СѓС‚ СЂР°Р·РЅС‹Рµ k).
    // -----------------------------------------------------------------------
    void color_parallel(const std::vector<size_t>& pair_i,
                        const std::vector<size_t>& pair_j,
                        size_t niac, size_t ntotal)
    {
        // Р’РµСЂС…РЅСЏСЏ РѕС†РµРЅРєР° С‡РёСЃР»Р° С†РІРµС‚РѕРІ: СЃС‚РµРїРµРЅСЊ РїР°СЂС‹ РІ РіСЂР°С„Рµ РєРѕРЅС„Р»РёРєС‚РѕРІ
        // <= (deg(i)-1)+(deg(j)-1) <= 2*(max_degree-1); Р¶Р°РґРЅР°СЏ РѕРєСЂР°СЃРєР° <= СЃС‚РµРїРµРЅСЊ+1.
        size_t max_degree = 0;
        for (size_t v = 0; v < ntotal; ++v)
        {
            size_t deg = pair_offset[v + 1] - pair_offset[v];
            if (deg > max_degree) max_degree = deg;
        }
        fb_stride = static_cast<int>(2 * max_degree + 1);
        if (fb_stride < 16) fb_stride = 16;  // РЅРёР¶РЅСЏСЏ РіСЂР°РЅРёС†Р° РѕС†РµРЅРєРё С†РІРµС‚РѕРІ

        int nthreads = omp_get_max_threads();
        if (static_cast<size_t>(fb_stride) * static_cast<size_t>(nthreads) > fb_par.size())
            fb_par.assign(static_cast<size_t>(nthreads) * static_cast<size_t>(fb_stride), -1);

        // NB: РїРѕСЂСЏРґРѕРє largest-degree-first Р·РґРµСЃСЊ РќР• РїСЂРёРјРµРЅСЏРµС‚СЃСЏ. РћРЅ СЂР°Р·СЂСѓС€Р°РµС‚
        // CSR-Р»РѕРєР°Р»СЊРЅРѕСЃС‚СЊ mark()/check() (С‡С‚РµРЅРёРµ pair_csr РїРѕРґ СЂР°Р·Р±СЂРѕСЃР°РЅРЅС‹РјРё РїРѕ
        // СЃС‚РµРїРµРЅСЏРј С‡Р°СЃС‚РёС†Р°РјРё РІРјРµСЃС‚Рѕ РїРѕС‚РѕРєРѕРІРѕРіРѕ СЃРєР°РЅРёСЂРѕРІР°РЅРёСЏ РїРѕРґ Р±Р»РѕРєРѕРј) Рё Р·Р°РјРµРґР»СЏРµС‚
        // С„Р°Р·С‹ (a)+(b) ~2Г—. fewer colors РґР»СЏ РїРѕС‚СЂРµР±РёС‚РµР»РµР№ РґР°С‘С‚ coalesce_colors,
        // РЅРµ Р·Р°РІРёСЃСЏС‰РµРµ РѕС‚ РїРѕСЂСЏРґРєР° РѕРєСЂР°СЃРєРё; largest-degree-first РѕСЃС‚Р°С‘С‚СЃСЏ РІ color_serial.

        double _a0 = omp_get_wtime();
        // ---- (a) РЎРїРµРєСѓР»СЏС‚РёРІРЅР°СЏ РїР°СЂР°Р»Р»РµР»СЊРЅР°СЏ РѕРєСЂР°СЃРєР° (Р»РѕРєР°Р»СЊРЅРѕ-РѕСЂРёРµРЅС‚РёСЂРѕРІР°РЅРЅР°СЏ) ----
        // РџРѕС‚РѕРєРё РґРµР»СЏС‚ С‡Р°СЃС‚РёС†С‹ РЅР° Р±Р»РѕРєРё [v_start, v_end); РєР°Р¶РґС‹Р№ РїРѕС‚РѕРє РєСЂР°СЃРёС‚ РїР°СЂС‹,
        // С‡СЊСЏ РџР•Р Р’РђРЇ С‡Р°СЃС‚РёС†Р° pair_i[k] Р»РµР¶РёС‚ РІ РµРіРѕ Р±Р»РѕРєРµ. РўРѕРіРґР° mark(pair_i[k])
        // С‡РёС‚Р°РµС‚ pair_csr РїРѕРґ С‡Р°СЃС‚РёС†Р°РјРё Р±Р»РѕРєР° вЂ” РїРѕС‚РѕРєРѕРІРѕРµ СЃРєР°РЅРёСЂРѕРІР°РЅРёРµ РЅРµРїСЂРµСЂС‹РІРЅРѕРіРѕ
        // СЃРµРіРјРµРЅС‚Р° pair_csr => РѕС‚Р»РёС‡РЅР°СЏ CSR-Р»РѕРєР°Р»СЊРЅРѕСЃС‚СЊ (~1.88Г— РїСЂРё P=8). РљР°Р¶РґР°СЏ
        // РїР°СЂР° РєСЂР°СЃРёС‚СЃСЏ СЂРѕРІРЅРѕ РѕРґРЅРёРј РїРѕС‚РѕРєРѕРј => Р·Р°РїРёСЃСЊ color[k] Р±РµСЃРєРѕРЅС„Р»РёРєС‚РЅР°.
        // fb_par вЂ” per-thread СЃСЂРµР·С‹ forbidden: fb_par[tid*fb_stride + c] == k.
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            int* fb = fb_par.data() + static_cast<size_t>(tid) * static_cast<size_t>(fb_stride);
            size_t v_start = (ntotal * static_cast<size_t>(tid)) / static_cast<size_t>(nthreads);
            size_t v_end = (ntotal * static_cast<size_t>(tid + 1)) / static_cast<size_t>(nthreads);
            auto mark = [&](size_t v, int cur_k) {
                for (size_t p = pair_offset[v]; p < pair_offset[v + 1]; ++p)
                {
                    size_t k2 = pair_csr[p];
                    if (k2 != static_cast<size_t>(cur_k))
                    {
                        int c2 = color[k2];  // РѕРїСЂРѕРјРµС‚С‡РёРІРѕРµ С‡С‚РµРЅРёРµ
                        if (c2 >= 0 && c2 < fb_stride)
                            fb[c2] = cur_k;
                    }
                }
            };
            for (size_t v = v_start; v < v_end; ++v)
            {
                for (size_t p = pair_offset[v]; p < pair_offset[v + 1]; ++p)
                {
                    size_t k = pair_csr[p];
                    if (pair_i[k] != v)
                        continue;  // РїР°СЂСѓ РєСЂР°СЃРёС‚ РїРѕС‚РѕРє Р±Р»РѕРєР° РµС‘ РїРµСЂРІРѕР№ С‡Р°СЃС‚РёС†С‹
                    size_t i = pair_i[k], j = pair_j[k];
                    const int marker = static_cast<int>(k);
                    mark(i, marker);
                    mark(j, marker);
                    int c = 0;
                    while (c < fb_stride && fb[c] == marker)
                        ++c;
                    color[k] = (c < fb_stride) ? c : -2;  // -2: РѕС†РµРЅРєР° РїСЂРµРІС‹С€РµРЅР° -> СЃРµСЂРёР№РЅРѕРµ СЂР°Р·СЂРµС€РµРЅРёРµ
                }
            }
        }
        stat_phase_a_us = static_cast<long>((omp_get_wtime() - _a0) * 1e6);
        double _b0 = omp_get_wtime();

        // ---- (b) РџР°СЂР°Р»Р»РµР»СЊРЅРѕРµ РѕР±РЅР°СЂСѓР¶РµРЅРёРµ РєРѕРЅС„Р»РёРєС‚РѕРІ (Р»РѕРєР°Р»СЊРЅРѕ-РѕСЂРёРµРЅС‚РёСЂРѕРІР°РЅРЅРѕРµ) ----
        // РћР±С…РѕРґ РїРѕ С‡Р°СЃС‚РёС†Р°Рј-Р±Р»РѕРєР°Рј (РєР°Рє С„Р°Р·Р° a): check(pair_i[k]) СЃРєР°РЅРёСЂСѓРµС‚ pair_csr
        // РїРѕРґ С‡Р°СЃС‚РёС†Р°РјРё Р±Р»РѕРєР° вЂ” РїРѕС‚РѕРєРѕРІРѕРµ С‡С‚РµРЅРёРµ (CSR-Р»РѕРєР°Р»СЊРЅРѕСЃС‚СЊ). РљР°Р¶РґСѓСЋ РїР°СЂСѓ
        // РїСЂРѕРІРµСЂСЏРµС‚ РїРѕС‚РѕРє Р±Р»РѕРєР° РµС‘ РїРµСЂРІРѕР№ С‡Р°СЃС‚РёС†С‹ (Р·Р°РїРёСЃСЊ conflict_flag Р±РµСЃРєРѕРЅС„Р»РёРєС‚РЅР°).
        // conflict_flag[k]=1, РµСЃР»Рё С†РІРµС‚ СЃРѕРІРїР°Р» СЃ СЃРѕСЃРµРґРѕРј Р»РёР±Рѕ == -2.
        std::fill_n(conflict_flag.begin(), niac, 0);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            size_t v_start = (ntotal * static_cast<size_t>(tid)) / static_cast<size_t>(nthreads);
            size_t v_end = (ntotal * static_cast<size_t>(tid + 1)) / static_cast<size_t>(nthreads);
            for (size_t v = v_start; v < v_end; ++v)
            {
                for (size_t p = pair_offset[v]; p < pair_offset[v + 1]; ++p)
                {
                    size_t k = pair_csr[p];
                    if (pair_i[k] != v)
                        continue;  // РїСЂРѕРІРµСЂСЏРµС‚ РїРѕС‚РѕРє Р±Р»РѕРєР° РїРµСЂРІРѕР№ С‡Р°СЃС‚РёС†С‹
                    int c = color[k];
                    if (c < 0)
                    {
                        conflict_flag[k] = 1;
                        continue;
                    }
                    size_t i = pair_i[k], j = pair_j[k];
                    bool conf = false;
                    auto check = [&](size_t vv) {
                        for (size_t pp = pair_offset[vv]; pp < pair_offset[vv + 1]; ++pp)
                        {
                            size_t k2 = pair_csr[pp];
                            if (k2 != k && color[k2] == c)
                            {
                                conf = true;
                                return;
                            }
                        }
                    };
                    check(i);
                    if (!conf) check(j);
                    if (conf) conflict_flag[k] = 1;
                }
            }
        }

        stat_phase_b_us = static_cast<long>((omp_get_wtime() - _b0) * 1e6);
        double _c0 = omp_get_wtime();

        // ---- (c) РЎРµСЂРёР№РЅРѕРµ СЂР°Р·СЂРµС€РµРЅРёРµ РєРѕРЅС„Р»РёРєС‚СѓСЋС‰РёС… РїР°СЂ ----
        // РџРѕРґСЃС‡С‘С‚ РєРѕРЅС„Р»РёРєС‚СѓСЋС‰РёС… РїР°СЂ (РґРёР°РіРЅРѕСЃС‚РёРєР°).
        stat_conflicts = 0;
        for (size_t _k = 0; _k < niac; ++_k)
            if (conflict_flag[_k]) ++stat_conflicts;
        // num_colors вЂ” РјР°РєСЃРёРјСѓРј С†РІРµС‚РѕРІ РёР· СЃРїРµРєСѓР»СЏС‚РёРІРЅРѕР№ С„Р°Р·С‹.
        num_colors = 0;
        for (size_t k = 0; k < niac; ++k)
            if (color[k] >= 0 && color[k] + 1 > num_colors)
                num_colors = color[k] + 1;

        std::fill(forbidden.begin(), forbidden.end(), -1);
        auto mark_neighbors = [&](size_t v, int cur_k) {
            for (size_t p = pair_offset[v]; p < pair_offset[v + 1]; ++p)
            {
                size_t k2 = pair_csr[p];
                if (k2 != static_cast<size_t>(cur_k) && color[k2] >= 0)
                    forbidden[static_cast<size_t>(color[k2])] = cur_k;
            }
        };

        // РЎРµСЂРёР№РЅС‹Р№ РїСЂРѕС…РѕРґ РїРѕ РІСЃРµРј РїР°СЂР°Рј (РєРѕРЅС„Р»РёРєС‚СѓСЋС‰РёС… РјР°Р»Рѕ РїСЂРё С‚РёРїРѕРІС‹С… SPH-РіСЂР°С„Р°С…).
        for (size_t k = 0; k < niac; ++k)
        {
            if (!conflict_flag[k] && color[k] >= 0)
                continue;  // РїР°СЂР° РєРѕСЂСЂРµРєС‚РЅР° вЂ” СЃРїРµРєСѓР»СЏС‚РёРІРЅС‹Р№ С†РІРµС‚ РѕСЃС‚Р°РІР»СЏРµРј

            size_t i = pair_i[k], j = pair_j[k];
            if (static_cast<int>(forbidden.size()) < num_colors + 1)
                forbidden.resize(static_cast<size_t>(num_colors) + 1, -1);

            const int marker = static_cast<int>(k);
            mark_neighbors(i, marker);
            mark_neighbors(j, marker);

            int c = 0;
            while (c < num_colors && forbidden[static_cast<size_t>(c)] == marker)
                ++c;
            if (c == num_colors)
            {
                num_colors++;
                forbidden.push_back(-1);
            }
            color[k] = c;
        }
        stat_phase_c_us = static_cast<long>((omp_get_wtime() - _c0) * 1e6);
    }

    // -----------------------------------------------------------------------
    // for_each_color вЂ” СѓРЅРёС„РёС†РёСЂРѕРІР°РЅРЅС‹Р№ РёРЅС‚РµСЂС„РµР№СЃ РїРѕС‚СЂРµР±РёС‚РµР»РµР№ РїР°СЂ.
    // Р”Р»СЏ РєР°Р¶РґРѕРіРѕ С†РІРµС‚Р° c РІС‹Р·С‹РІР°РµС‚ body(c, begin, end), РіРґРµ [begin, end) вЂ”
    // РґРёР°РїР°Р·РѕРЅ РёРЅРґРµРєСЃРѕРІ РІ pairs_by_color, СЃРѕРґРµСЂР¶Р°С‰РёР№ РїР°СЂС‹ С†РІРµС‚Р° c.
    // Р’РЅСѓС‚СЂРё body РїРѕС‚СЂРµР±РёС‚РµР»СЊ Р·Р°РїСѓСЃРєР°РµС‚ С‡РёСЃС‚С‹Р№ #pragma omp parallel for
    // РїРѕ РґРёР°РїР°Р·РѕРЅСѓ [begin, end) Р‘Р•Р— atomic (РїР°СЂС‹ РѕРґРЅРѕРіРѕ С†РІРµС‚Р° РЅРµ РїРµСЂРµСЃРµРєР°СЋС‚СЃСЏ).
    // РњРµР¶РґСѓ С†РІРµС‚Р°РјРё OpenMP-Р±Р°СЂСЊРµСЂ РѕР±СЏР·Р°С‚РµР»РµРЅ (С†РІРµС‚ c+1 РјРѕР¶РµС‚ СЂР°Р·РґРµР»СЏС‚СЊ С‡Р°СЃС‚РёС†Сѓ
    // СЃ С†РІРµС‚РѕРј c).
    // -----------------------------------------------------------------------
    template <class Body>
    void for_each_color(Body&& body) const
    {
        for (int c = 0; c < num_colors; ++c)
            body(c, color_offset[c], color_offset[c + 1]);
    }

    // -----------------------------------------------------------------------
    // is_valid_coloring вЂ” РїСЂРѕРІРµСЂРєР° РёРЅРІР°СЂРёР°РЅС‚Р° РїСЂР°РІРёР»СЊРЅРѕР№ СЂР°СЃРєСЂР°СЃРєРё (РґР»СЏ С‚РµСЃС‚РѕРІ).
    // Р’РѕР·РІСЂР°С‰Р°РµС‚ true, РµСЃР»Рё РЅРёРєР°РєРёРµ РґРІРµ РїР°СЂС‹ РѕРґРЅРѕРіРѕ С†РІРµС‚Р° РЅРµ СЂР°Р·РґРµР»СЏСЋС‚ С‡Р°СЃС‚РёС†Сѓ.
    // -----------------------------------------------------------------------
    bool is_valid_coloring(const std::vector<size_t>& pair_i,
                           const std::vector<size_t>& pair_j,
                           size_t niac) const
    {
        if (num_colors <= 0)
            return niac == 0;
        for (int c = 0; c < num_colors; ++c)
        {
            // РЎРѕР±СЂР°С‚СЊ РІСЃРµ С‡Р°СЃС‚РёС†С‹, Р·Р°С‚СЂРѕРЅСѓС‚С‹Рµ РїР°СЂР°РјРё С†РІРµС‚Р° c, Рё РїСЂРѕРІРµСЂРёС‚СЊ СѓРЅРёРєР°Р»СЊРЅРѕСЃС‚СЊ.
            std::vector<size_t> touched;
            for (size_t p = color_offset[static_cast<size_t>(c)]; p < color_offset[static_cast<size_t>(c) + 1]; ++p)
            {
                size_t k = pairs_by_color[p];
                touched.push_back(pair_i[k]);
                touched.push_back(pair_j[k]);
            }
            std::sort(touched.begin(), touched.end());
            for (size_t t = 1; t < touched.size(); ++t)
            {
                if (touched[t] == touched[t - 1])
                    return false;  // РѕРґРЅР° С‡Р°СЃС‚РёС†Р° РІ РґРІСѓС… РїР°СЂР°С… С†РІРµС‚Р° c
            }
        }
        return true;
    }
};
