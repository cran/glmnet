#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/binomial_base.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::binomial, 
    util::mode_type<util::glm_type::binomial>::two_class,
    ElnetPointInternalPolicy>
        : ElnetPointBinomialBase<
            ElnetPoint<
                util::glm_type::binomial, 
                util::mode_type<util::glm_type::binomial>::two_class,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointBinomialBase<
            ElnetPoint<
                util::glm_type::binomial, 
                util::mode_type<util::glm_type::binomial>::two_class,
                ElnetPointInternalPolicy> >;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
