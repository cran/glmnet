#pragma once
#include <glmnetpp_bits/elnet_point/gaussian_wls.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::gaussian, 
    util::mode_type<util::glm_type::gaussian>::wls,
    ElnetPointInternalPolicy>
        : ElnetPoint<
                util::glm_type::gaussian, 
                util::mode_type<util::glm_type::gaussian>::wls,
                ElnetPointInternalPolicy>
{
private:
    using base_t = ElnetPoint<util::glm_type::gaussian, 
                   util::mode_type<util::glm_type::gaussian>::wls,
                   ElnetPointInternalPolicy>;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
