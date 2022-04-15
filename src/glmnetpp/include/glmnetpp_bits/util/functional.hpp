#pragma once

namespace glmnetpp {
namespace util {
namespace details {

template <bool cond>
struct IfElseFunc
{
    template <class F1, class F2>
    constexpr inline static auto eval(F1 f1, F2)
    {
        return f1();
    }
};

template <>
struct IfElseFunc<false>
{
    template <class F1, class F2>
    constexpr inline static auto eval(F1, F2 f2)
    {
        return f2();
    }
};

} // namespace details

/*
 * Compile-time if-else based on condition cond.
 * If cond is true, f1 is executed.
 * Otherwise, f2 is executed.
 */
template <bool cond, class F1, class F2>
constexpr inline auto if_else(F1 f1, F2 f2)
{
    return details::IfElseFunc<cond>::eval(f1, f2);
}

} // namespace util
} // namespace glmnetpp
