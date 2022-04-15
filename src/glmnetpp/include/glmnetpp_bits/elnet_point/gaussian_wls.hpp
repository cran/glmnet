#pragma once
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_base.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::gaussian, 
    util::mode_type<util::glm_type::gaussian>::wls,
    ElnetPointInternalPolicy>
        : ElnetPointGaussianBase<
            ElnetPoint<
                util::glm_type::gaussian, 
                util::mode_type<util::glm_type::gaussian>::wls,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointGaussianBase<
        ElnetPoint<util::glm_type::gaussian, 
                   util::mode_type<util::glm_type::gaussian>::wls,
                   ElnetPointInternalPolicy> >;
    using typename base_t::update_t;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    
    // hide some unnecessary parts of the base class
    using base_t::fit;

    // dummy struct to stay consistent with interface.
    struct PointConfigPack {};

public:
    using base_t::base_t;

    template <class IntType>
    GLMNETPP_STRONG_INLINE
    void fit(IntType m, IntType& jerr) 
    {
        try {
            base_t::fit(PointConfigPack());
        } 
        catch (const util::elnet_error& e) {
            jerr = e.err_code(m-1); // m is 1-indexed
        }
    }

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
