// test_fiber_dim_guard.cpp
// Проверяет ветку "dim != 3" границы Р9 (notes/fiber-model/fiber-plan.md, этап Ф2, пункт 5):
// fiber::check_fiber_requires_dim3() обязана бросать std::runtime_error,
// когда сборка сделана не с SPH_DIM=3.
//
// dim — compile-time constexpr (param.hpp, ADR-002). Эта цель СОЗНАТЕЛЬНО
// собирается БЕЗ SPH_DIM (см. CMakeLists.txt): dim = 1 по умолчанию, поэтому
// ветка if constexpr (dim != 3) внутри check_fiber_requires_dim3() реально
// компилируется и исполняется здесь. В test_fiber_model (SPH_DIM=3) эта же
// функция вызывается на непротиворечивость (не бросает при dim == 3), но
// ветка throw там недостижима в принципе — см. комментарий в
// fiber/fiber_system.hpp.
#include "fiber/fiber_system.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>

static_assert(dim != 3, "test_fiber_dim_guard должен собираться БЕЗ SPH_DIM=3 (проверяет именно ветку dim != 3)");

int main()
{
    int failures = 0;

    bool threw = false;
    std::string message;
    try
    {
        fiber::check_fiber_requires_dim3();
    }
    catch (const std::runtime_error& e)
    {
        threw = true;
        message = e.what();
    }

    if (!threw)
    {
        std::printf("FAIL: check_fiber_requires_dim3() не бросила исключение при dim=%zu\n", dim);
        ++failures;
    }
    else
    {
        std::printf("ok: check_fiber_requires_dim3() бросила: %s\n", message.c_str());
        if (message.find("SPH_DIM=3") == std::string::npos && message.find("dim") == std::string::npos)
        {
            std::printf("FAIL: сообщение об ошибке не упоминает dim/SPH_DIM: %s\n", message.c_str());
            ++failures;
        }
    }

    if (failures == 0)
    {
        std::printf("test_fiber_dim_guard: OK\n");
    }
    else
    {
        std::printf("test_fiber_dim_guard: %d failure(s)\n", failures);
    }
    return failures;
}
