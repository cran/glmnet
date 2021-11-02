#pragma once
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/type_traits.hpp>

namespace glmnetpp {
namespace details {

template <class T> struct traits;

template <util::glm_type g
        , util::mode_type<g> m
        , class ElnetPointInternalType>
struct traits<ElnetPoint<g, m, ElnetPointInternalType> >
{
    static constexpr util::glm_type glm = g;
    static constexpr util::mode_type<g> mode = m;
    using internal_t = ElnetPointInternalType;
};

template <util::glm_type g
        , util::mode_type<g> m
        , class ElnetPointInternalType>
struct traits<SpElnetPoint<g, m, ElnetPointInternalType> >
{
    static constexpr util::glm_type glm = g;
    static constexpr util::mode_type<g> mode = m;
    using internal_t = ElnetPointInternalType;
};

} // namespace details
} // namespace glmnetpp
