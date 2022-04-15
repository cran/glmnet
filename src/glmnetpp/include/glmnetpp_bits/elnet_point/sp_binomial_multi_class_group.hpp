#pragma once
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/binomial_multi_class_group.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::binomial, 
    util::mode_type<util::glm_type::binomial>::multi_class_group,
    ElnetPointInternalPolicy>
        : ElnetPoint<
                util::glm_type::binomial, 
                util::mode_type<util::glm_type::binomial>::multi_class_group,
                ElnetPointInternalPolicy>
{
private:
    using base_t = ElnetPoint<
        util::glm_type::binomial, 
        util::mode_type<util::glm_type::binomial>::multi_class_group,
        ElnetPointInternalPolicy>;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
