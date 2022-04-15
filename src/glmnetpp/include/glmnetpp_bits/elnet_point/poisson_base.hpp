#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/base.hpp>

namespace glmnetpp {

template <class Derived>
struct ElnetPointPoissonBase : 
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
