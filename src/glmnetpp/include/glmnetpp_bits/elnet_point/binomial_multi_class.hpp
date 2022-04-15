#pragma once
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/binomial_base.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::binomial, 
    util::mode_type<util::glm_type::binomial>::multi_class,
    ElnetPointInternalPolicy>
        : ElnetPointBinomialBase<
            ElnetPoint<
                util::glm_type::binomial, 
                util::mode_type<util::glm_type::binomial>::multi_class,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointBinomialBase<
        ElnetPoint<util::glm_type::binomial, 
                   util::mode_type<util::glm_type::binomial>::multi_class,
                   ElnetPointInternalPolicy> >;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    void irls(const PointConfigPack& pack) 
    {
        this->initialize(pack);

        while (1) {
            this->reset_converged();
            if (this->has_reached_max_passes()) {
                throw util::maxit_reached_error();
            }
            try {
                std::for_each(
                        this->class_begin(),
                        this->class_end(),
                        [&, this](auto ic) {
                            state_t state = this->setup_wls(ic);
                            if (state == state_t::continue_) return;
                            base_t::wls(pack);
                            this->update_irls_class();
                        });
            }
            // catch so that we can still finish some of the IRLS update
            catch (const util::max_active_reached_error& e) {}
                
            state_t state = this->update_irls(pack);
            if (state == state_t::break_) break;
        }
    }
};

} // namespace glmnetpp
