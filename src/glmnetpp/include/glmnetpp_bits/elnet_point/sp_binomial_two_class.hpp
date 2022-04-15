#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/binomial_base.hpp>

namespace glmnetpp {

template <class SpElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::binomial, 
    util::mode_type<util::glm_type::binomial>::two_class,
    SpElnetPointInternalPolicy>
        : ElnetPoint<
            util::glm_type::binomial, 
            util::mode_type<util::glm_type::binomial>::two_class,
            SpElnetPointInternalPolicy>
{
private:
    using base_t = ElnetPoint<
                util::glm_type::binomial, 
                util::mode_type<util::glm_type::binomial>::two_class,
                SpElnetPointInternalPolicy>;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
