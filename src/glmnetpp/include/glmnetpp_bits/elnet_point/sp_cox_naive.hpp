#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/base.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_cox_naive.hpp>

namespace glmnetpp {

/*
 * Sparse Cox point solver (naive method).
 * This is a direct implementation inheriting from the CRTP base,
 * bypassing ElnetPoint<cox> to avoid constructor inheritance issues
 * with the different parameter count (16 vs 14).
 */
template <class SpElnetPointInternalPolicy>
struct SpElnetPoint<
    util::glm_type::cox,
    util::mode_type<util::glm_type::cox>::naive,
    SpElnetPointInternalPolicy>
        : ElnetPointNonLinearCRTPBase<
            SpElnetPoint<
                util::glm_type::cox,
                util::mode_type<util::glm_type::cox>::naive,
                SpElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointNonLinearCRTPBase<
            SpElnetPoint<
                util::glm_type::cox,
                util::mode_type<util::glm_type::cox>::naive,
                SpElnetPointInternalPolicy> >;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;
};

} // namespace glmnetpp
