#pragma once
#include <glmnetpp_bits/util/types.hpp>

namespace glmnetpp {

template <util::glm_type g
        , util::mode_type<g> mode
        , class ValueType = double
        , class IndexType = int
        , class BoolType = bool>
struct ElnetPointInternal;

template <util::glm_type g
        , util::mode_type<g> mode
        , class ValueType = double
        , class IndexType = int
        , class BoolType = bool>
struct SpElnetPointInternal;

} // namespace glmnetpp
