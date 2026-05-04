#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/base.hpp>

namespace glmnetpp {

/*
 * CRTP base class for Cox point solvers.
 * Cox uses IRLS like Poisson, but with no intercept.
 */
template <class Derived>
struct ElnetPointCoxBase :
    ElnetPointNonLinearCRTPBase<Derived>
{
private:
    using base_t = ElnetPointNonLinearCRTPBase<Derived>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
