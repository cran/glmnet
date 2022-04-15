#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/poisson_base.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::poisson, 
    util::mode_type<util::glm_type::poisson>::naive,
    ElnetPointInternalPolicy>
        : ElnetPointPoissonBase<
            ElnetPoint<
                util::glm_type::poisson, 
                util::mode_type<util::glm_type::poisson>::naive,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointPoissonBase<
            ElnetPoint<
                util::glm_type::poisson, 
                util::mode_type<util::glm_type::poisson>::naive,
                ElnetPointInternalPolicy> >;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
