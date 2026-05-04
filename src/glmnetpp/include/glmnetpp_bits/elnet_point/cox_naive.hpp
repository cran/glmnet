#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/cox_base.hpp>

namespace glmnetpp {

/*
 * Cox point solver (naive/dense method).
 * Specialization that delegates to the internal implementation.
 */
template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::cox,
    util::mode_type<util::glm_type::cox>::naive,
    ElnetPointInternalPolicy>
        : ElnetPointCoxBase<
            ElnetPoint<
                util::glm_type::cox,
                util::mode_type<util::glm_type::cox>::naive,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointCoxBase<
            ElnetPoint<
                util::glm_type::cox,
                util::mode_type<util::glm_type::cox>::naive,
                ElnetPointInternalPolicy> >;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
