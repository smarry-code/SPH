#ifndef RHS_WORKSPACE_HPP
#define RHS_WORKSPACE_HPP

#include <any>
#include <string>
#include <unordered_map>

// ============================================================================
// rhs/rhs_workspace.hpp — RhsWorkspace: контейнер временной памяти ОДНОГО
// вычисления RHS (notes/fiber-model/fiber-problem-statement.md, раздел "Временная память"; notes/fiber-model/fiber-plan.md
// §3.5): пространственные индексы, буферы пар, primitives, contributions,
// шарды reduce, производные величины Solid, виртуальные частицы.
//
// "Память переиспользуется между вызовами RHS и только растёт" (fiber-plan.md):
// RhsWorkspace никогда сам не освобождает уже выделенную память — вызывающий
// код (например, ContributionBuffer::clear() перед новым RHS) сам решает,
// когда логически "обнулить" содержимое буфера, сохранив его capacity.
//
// На этапе Ф1 конкретные типы (SolidDelta, FiberState, GhostSet,
// SpatialIndex<Fiber-конкретной-размерности> и т.п.) ещё не существуют —
// задача прямо запрещает физику. Поэтому RhsWorkspace реализован как
// типобезопасный реестр именованных переиспользуемых буферов произвольного
// типа: каждый следующий этап (Ф2+) получает свой переиспользуемый буфер
// через get<T>("имя"), не требуя изменений здесь же в RhsWorkspace.
// ============================================================================

namespace rhs
{

class RhsWorkspace
{
public:
    // Возвращает переиспользуемый буфер типа T с данным ключом. При первом
    // обращении создаётся T{} и сохраняется; при повторных обращениях (в том
    // числе из следующего вызова RHS) возвращается тот же объект — тот же
    // адрес, та же уже выделенная память.
    template <class T>
    T& get(const std::string& key)
    {
        auto it = storage_.find(key);
        if (it == storage_.end())
        {
            it = storage_.emplace(key, std::any(T{})).first;
        }
        return std::any_cast<T&>(it->second);
    }

    bool contains(const std::string& key) const { return storage_.find(key) != storage_.end(); }

    // Полный сброс workspace (освобождает всю память). Это НЕ вызывается
    // между последовательными вычислениями RHS одного расчёта (там память
    // обязана только расти) — предназначено для независимых прогонов/тестов.
    void reset() { storage_.clear(); }

private:
    std::unordered_map<std::string, std::any> storage_;
};

} // namespace rhs

#endif // RHS_WORKSPACE_HPP
