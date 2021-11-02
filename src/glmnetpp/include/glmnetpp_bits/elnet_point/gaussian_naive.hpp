#pragma once
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_base.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::gaussian, 
    util::mode_type<util::glm_type::gaussian>::naive,
    ElnetPointInternalPolicy>
        : ElnetPointGaussianBase<
            ElnetPoint<
                util::glm_type::gaussian, 
                util::mode_type<util::glm_type::gaussian>::naive,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointGaussianBase<
        ElnetPoint<util::glm_type::gaussian, 
                   util::mode_type<util::glm_type::gaussian>::naive,
                   ElnetPointInternalPolicy> >;
    using typename base_t::update_type;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using base_t::partial_fit;
    using base_t::update;

public:
    using base_t::base_t;
    using base_t::abs_grad;

    void partial_fit(index_t m, value_t ab, value_t dem)
    {
        base_t::partial_fit(m, ab, dem);

        // For some reason, naive version does this extra check.
        if (this->has_reached_max_passes()) { 
            throw util::maxit_reached_error();
        }
    }

    template <update_type upd>
    void update(index_t k, value_t ab, value_t dem)
    {
        value_t beta_diff = 0;

        auto state = base_t::template update<upd>(k, ab, dem, beta_diff);
        if (state == state_t::continue_) return;
        
        this->update_resid(k, beta_diff);
    }
};

} // namespace glmnetpp
