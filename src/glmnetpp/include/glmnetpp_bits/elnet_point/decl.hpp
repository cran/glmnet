#pragma once
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>

namespace glmnetpp {

template <class ElnetPointDerived>
struct ElnetPointCRTPBase;
template <class ElnetPointDerived>
struct ElnetPointNonLinearCRTPBase;

template <class Derived>
struct ElnetPointBinomialBase;
template <class Derived>
struct ElnetPointGaussianBase;
template <class Derived>
struct ElnetPointPoissonBase; 

template <util::glm_type glm
        , util::mode_type<glm> mode
        , class ElnetPointInternalPolicy=ElnetPointInternal<glm, mode> >
struct ElnetPoint;

template <util::glm_type glm
        , util::mode_type<glm> mode
        , class SpElnetPointInternalPolicy=SpElnetPointInternal<glm, mode> >
struct SpElnetPoint;


} // namespace glmnetpp
