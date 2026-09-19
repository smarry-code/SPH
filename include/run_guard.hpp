#pragma once
// run_guard.hpp
// Сторож расчёта: обрывает главный цикл по времени, если расчёт исчерпал
// отведённое время или перестал продвигаться (коллапс шага).
//
// Зачем. time_integration() — монолитный цикл `while (time < t_total)`. Если шаг
// по времени схлопывается (частица с аномальной скоростью звука или нулевой
// длиной сглаживания), цикл не зависает в смысле ОС: он продолжает считать, но
// модельное время практически стоит. Снаружи это неотличимо от долгого расчёта,
// и процесс живёт до тех пор, пока его не снимут вручную.
//
// Сторож различает ДВА исхода, потому что чинить их надо по-разному:
//   WallClock — расчёт шёл штатно, шаг здоров, просто не уложился в бюджет.
//               Лечится увеличением бюджета или огрублением сетки.
//   Stalled   — шаг деградировал настолько, что расчёт не финиширует никогда.
//               Это дефект модели или данных, бюджет тут не поможет.
//
// По умолчанию выключен (лимиты нулевые): поведение конвейера без явной
// настройки не меняется, main.cpp сторож не передаёт.
#include <chrono>
#include <cstddef>
#include <sstream>
#include <string>

class RunGuard
{
public:
    enum class Reason
    {
        None,       // не сработал
        WallClock,  // исчерпан бюджет времени
        Stalled     // расчёт перестал продвигаться
    };

    RunGuard() = default;

    // wall_limit_s        — бюджет времени на весь расчёт, с; <= 0 — без лимита.
    // max_projected_steps — потолок оценки «сколько шагов осталось при текущем
    //                       dt»; превышение, продержавшееся patience шагов
    //                       подряд, считается застоем. 0 — не следить.
    // patience            — сколько шагов подряд терпеть превышение.
    RunGuard(double wall_limit_s, double max_projected_steps, size_t patience)
        : wall_limit_s_(wall_limit_s),
          max_projected_steps_(max_projected_steps),
          patience_(patience)
    {
    }

    // Вызывается ровно один раз за шаг, после обновления time.
    // Возвращает true, если цикл нужно прервать.
    bool should_stop(double dt, double time, double t_total)
    {
        ++steps_;
        if (dt_first_ <= 0.0 && dt > 0.0)
        {
            dt_first_ = dt;
        }
        if (dt > 0.0 && (dt_min_ <= 0.0 || dt < dt_min_))
        {
            dt_min_ = dt;
        }
        dt_last_  = dt;
        time_     = time;
        t_total_  = t_total;

        // Оценка остатка при текущем шаге. dt <= 0 — вырожденный случай:
        // время не растёт вовсе, остаток бесконечен.
        const double remaining = t_total - time;
        projected_ = (dt > 0.0) ? (remaining > 0.0 ? remaining / dt : 0.0)
                                : kInfProjection;

        if (max_projected_steps_ > 0.0 && projected_ > max_projected_steps_)
        {
            ++over_;
            if (over_ >= patience_)
            {
                reason_ = Reason::Stalled;
                return true;
            }
        }
        else
        {
            over_ = 0;
        }

        if (wall_limit_s_ > 0.0 && elapsed_s() > wall_limit_s_)
        {
            reason_ = Reason::WallClock;
            return true;
        }
        return false;
    }

    bool   triggered() const { return reason_ != Reason::None; }
    Reason reason()    const { return reason_; }

    double elapsed_s() const
    {
        return std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - started_)
            .count();
    }

    double progress()   const { return t_total_ > 0.0 ? time_ / t_total_ : 0.0; }
    double dt_first()   const { return dt_first_; }
    double dt_last()    const { return dt_last_; }
    double dt_min()     const { return dt_min_; }
    double projected()  const { return projected_; }
    size_t steps()      const { return steps_; }

    // Во сколько раз шаг просел относительно первого. 0 — оценить нечем.
    double dt_drop() const
    {
        return (dt_first_ > 0.0 && dt_last_ > 0.0) ? dt_first_ / dt_last_ : 0.0;
    }

    // Признак здорового шага на момент обрыва: расчёт реально двигался.
    // Порог 10x — шаг в ударной задаче законно проседает в разы, но не на
    // порядок относительно стартового.
    bool step_healthy() const
    {
        const double d = dt_drop();
        return d > 0.0 && d < 10.0;
    }

    // Готовая диагностика: что произошло и что с этим делать.
    std::string message() const
    {
        std::ostringstream os;
        os.setf(std::ios::scientific);
        os.precision(3);
        if (reason_ == Reason::WallClock)
        {
            os << "РАСЧЁТ ПРЕРВАН: исчерпан бюджет времени ("
               << fixed2(wall_limit_s_) << " с). ";
            if (step_healthy())
            {
                os << "Шаг по времени здоров (просел в " << fixed2(dt_drop())
                   << "x от стартового) — расчёт шёл штатно, но не уложился. "
                   << "Увеличьте бюджет или огрубите сетку.";
            }
            else
            {
                os << "Шаг по времени деградировал (просел в " << fixed2(dt_drop())
                   << "x от стартового) — вероятен коллапс шага, бюджет не поможет.";
            }
        }
        else if (reason_ == Reason::Stalled)
        {
            os << "РАСЧЁТ ПРЕРВАН: расчёт перестал продвигаться. "
               << "При текущем шаге до финиша осталось бы ~"
               << projected_ << " шагов (потолок " << max_projected_steps_
               << "). Это коллапс шага, а не нехватка времени.";
        }
        else
        {
            return "сторож не срабатывал";
        }
        os << "\n  пройдено " << fixed2(progress() * 100.0) << "% времени задачи (t = "
           << time_ << " из " << t_total_ << " с) за " << steps_ << " шагов и "
           << fixed2(elapsed_s()) << " с;"
           << "\n  dt: стартовый " << dt_first_ << ", последний " << dt_last_
           << ", минимальный " << dt_min_ << ".";
        return os.str();
    }

private:
    static std::string fixed2(double v)
    {
        std::ostringstream os;
        os.setf(std::ios::fixed);
        os.precision(2);
        os << v;
        return os.str();
    }

    static constexpr double kInfProjection = 1e300;

    double wall_limit_s_        = 0.0;
    double max_projected_steps_ = 0.0;
    size_t patience_            = 1;

    std::chrono::steady_clock::time_point started_ = std::chrono::steady_clock::now();
    Reason reason_   = Reason::None;
    double dt_first_ = 0.0;
    double dt_last_  = 0.0;
    double dt_min_   = 0.0;
    double time_     = 0.0;
    double t_total_  = 0.0;
    double projected_ = 0.0;
    size_t steps_    = 0;
    size_t over_     = 0;
};
