// checkpoint.hpp
// Бинарный чекпоинт для побитового рестарта расчёта (src/checkpoint.cpp).
//
// Сохраняется динамическое состояние Model (Solid на [0,ntotal), нити, топология
// связей с BondState, зажимы), состояние цикла и история решателя, а также
// SimulationConfig. Сцену (размеры, материалы нитей) вызывающий строит заново
// из аргументов запуска; read() проверяет, что она совпадает по размерам,
// и перезаписывает динамические поля.
//
// Формат: little-endian, magic "SPHCKPT\0", версия формата, SPH_DIM, флаги
// сборки, sizeof(double), блоки с тегом и длиной, в конце CRC32 всего файла.
#pragma once
#include <cstddef>
#include <string>
#include <vector>
#include "fiber/model.hpp"
#include "param.hpp"
#include "simulation_config.hpp"

namespace checkpoint
{

// Состояние цикла run_time_integration и RunGuard.
struct LoopState
{
    double time = 0.0, dt = 0.0, next_output_time = 0.0;
    size_t step = 0;                 // current_ts
    size_t guard_steps = 0;          // RunGuard::steps_
    double guard_dt_first = 0.0, guard_dt_min = 0.0;
};

// Состояние EnergyLogger; в файл пишется всегда, используется только при SPH_ENERGY_LOG.
struct EnergyState
{
    double chan[12] = {};
    double e_mech0 = 0.0;
    bool have_e0 = false;
};

// Копии истории PredictorCorrectorSolver (плоские векторы; пустой вектор допустим,
// например для EulerSolver, иначе длина ntotal / nfiber).
struct SolverHistory
{
    bool first_step = true;
    std::vector<double> du, drho;
    std::vector<double> dsxx, dsyy, dszz, dsxy, dsyz, dsxz;
    std::vector<double> dvx[dim];      // solid_rhs.dvx[d]
    std::vector<double> fiber_dv[dim]; // fiber_rhs.dv[d]
};

struct RestoreOptions
{
    bool force = false; // --force: разрешить менять настройки физики
};

// Различия сохранённых и новых настроек, строки "поле: старое -> новое".
struct ConfigDiff
{
    std::vector<std::string> allowed;   // параметры запуска (t_total, save_t, ...)
    std::vector<std::string> forbidden; // физика/численные настройки
};

struct Payload
{
    LoopState loop;
    SolverHistory hist;
    EnergyState energy;
    bool has_energy = false;
};

// Атомарная запись: временный файл + rename; существующий файл перед этим
// переименовывается в path + ".prev".
void write(const std::string& path, const SimulationConfig& cfg, const fiber::Model& model, size_t ntotal,
           const Payload& payload);

// Чтение. cfg: на входе настройки нового запуска, на выходе действующие
// (всегда новые значения; сохранённые нужны только для сравнения). model должна
// быть уже построена из аргументов запуска: размеры проверяются, динамические
// поля перезаписываются. ntotal на входе -- ожидаемое число частиц Solid.
// Бросает std::runtime_error (по-русски) при плохом magic/версии/SPH_DIM/флагах
// сборки/размерах/CRC/обрезанном файле и при запрещённом отличии настроек без force.
void read(const std::string& path, SimulationConfig& cfg, fiber::Model& model, size_t& ntotal, Payload& payload,
          const RestoreOptions& opt, ConfigDiff* diff = nullptr);

// Описание рисков --force (многострочный текст).
std::string force_risks_text();

} // namespace checkpoint
