#pragma once
#include <glmnetpp_bits/elnet_point/decl.hpp>

namespace glmnetpp {

template <class SpElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::gaussian, 
    util::mode_type<util::glm_type::gaussian>::naive,
    SpElnetPointInternalPolicy>
        : ElnetPoint<
            util::glm_type::gaussian, 
            util::mode_type<util::glm_type::gaussian>::naive,
            SpElnetPointInternalPolicy>
{
private:
    using base_t = ElnetPoint<
            util::glm_type::gaussian, 
            util::mode_type<util::glm_type::gaussian>::naive,
            SpElnetPointInternalPolicy>;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
