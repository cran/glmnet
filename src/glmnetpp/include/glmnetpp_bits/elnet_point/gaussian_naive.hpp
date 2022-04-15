#pragma once
#include <glmnetpp_bits/util/macros.hpp>
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
    using typename base_t::update_t;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;

public:
    using base_t::base_t;

    template <update_t upd, class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update(index_t k, const PointPackType& pack)
    {
        value_t beta_diff = 0;

        auto state = base_t::template update<upd>(k, pack, beta_diff);
        if (state == state_t::continue_) return;
        
        this->update_resid(k, beta_diff);
    }
};

} // namespace glmnetpp
