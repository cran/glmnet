#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/poisson_base.hpp>

namespace glmnetpp {

template <class SpElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::poisson, 
    util::mode_type<util::glm_type::poisson>::naive,
    SpElnetPointInternalPolicy>
        : ElnetPoint<
            util::glm_type::poisson, 
            util::mode_type<util::glm_type::poisson>::naive,
            SpElnetPointInternalPolicy>
{
private:
    using base_t = ElnetPoint<
                util::glm_type::poisson, 
                util::mode_type<util::glm_type::poisson>::naive,
                SpElnetPointInternalPolicy>;
public:
    using base_t::base_t;
};

} // namespace glmnetpp
