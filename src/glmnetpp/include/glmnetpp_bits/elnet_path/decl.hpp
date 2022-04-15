#pragma once
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>

namespace glmnetpp {

template <class ElnetPathDerived>
struct ElnetPathCRTPBase;

struct ElnetPathBase;
struct ElnetPathGaussianBase;
struct ElnetPathBinomialBase;

template <util::glm_type glm
        , util::mode_type<glm> mode
        , class ElnetPointPolicy=ElnetPoint<glm, mode> >
struct ElnetPath;

template <util::glm_type glm
        , util::mode_type<glm> mode
        , class SpElnetPointPolicy=SpElnetPoint<glm, mode> >
struct SpElnetPath;

} // namespace glmnetpp
