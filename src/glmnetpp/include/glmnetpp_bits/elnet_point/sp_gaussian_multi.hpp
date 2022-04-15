#pragma once
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_multi.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::gaussian, 
    util::mode_type<util::glm_type::gaussian>::multi,
    ElnetPointInternalPolicy>
        : ElnetPoint<
                util::glm_type::gaussian, 
                util::mode_type<util::glm_type::gaussian>::multi,
                ElnetPointInternalPolicy> 
{
private:
    using base_t = ElnetPoint<
                util::glm_type::gaussian, 
                util::mode_type<util::glm_type::gaussian>::multi,
                ElnetPointInternalPolicy>;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
