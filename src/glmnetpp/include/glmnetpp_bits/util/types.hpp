#pragma once

namespace glmnetpp {
namespace util {

// Types of glm supported.
enum class glm_type
{
    gaussian,
    binomial,
    poisson
};


// A struct that defines the enum class of mode types
// for a specific glm.
template <glm_type glm>
struct Mode;

// Specializations
template <>
struct Mode<glm_type::gaussian>
{
    enum class type
    {
        naive,
        cov,
        multi,
        wls
    };
};

template <>
struct Mode<glm_type::binomial>
{
    enum class type
    {
        two_class,
        multi_class,
        multi_class_group
    };
};

template <>
struct Mode<glm_type::poisson>
{
    enum class type
    {
        naive
    };
};

// Helper alias to get the mode enum class type for each glm.
template <glm_type glm>
using mode_type = typename Mode<glm>::type;

// Represents a state of a function as a way for the caller to do control flow.
enum class control_flow
{
    noop_,
    continue_,
    break_,
    return_
};

enum class update_type
{
    full,
    partial
};

} // namespace util
} // namespace glmnetpp
