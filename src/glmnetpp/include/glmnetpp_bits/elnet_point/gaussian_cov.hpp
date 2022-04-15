#pragma once
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_base.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/functional.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::gaussian, 
    util::mode_type<util::glm_type::gaussian>::cov,
    ElnetPointInternalPolicy>
        : ElnetPointGaussianBase<
            ElnetPoint<
                util::glm_type::gaussian, 
                util::mode_type<util::glm_type::gaussian>::cov,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointGaussianBase<
        ElnetPoint<util::glm_type::gaussian, 
                   util::mode_type<util::glm_type::gaussian>::cov,
                   ElnetPointInternalPolicy> >;
    using typename base_t::update_t;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;

public:
    using base_t::base_t;

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void partial_fit(const PointPackType& pack)
    {
        this->compress_active();

        base_t::partial_fit(pack);

        // update gradient to leave invariant
        this->update_compressed_active();
        this->update_grad_compressed_active();
    }

    template <update_t upd, class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update(index_t k, const PointPackType& pack)
    {
        value_t beta_diff = 0;

        auto state = base_t::template update<upd>(k, pack, beta_diff);
        if (state == state_t::continue_) return;

        // update gradient
        util::if_else<upd == update_t::full>(
                [=]() {
                    std::for_each(
                            this->all_begin(),
                            this->all_end(),
                            [=](auto j) {
                                if (!this->is_excluded(j)) {
                                    this->update_grad(j, k, beta_diff);
                                }
                            });
                },
                [=]() {
                    std::for_each(
                            this->active_begin(),
                            this->active_end(),
                            [=](auto j) {
                                this->update_grad(j, k, beta_diff);
                            });
                });
    }
};

} // namespace glmnetpp
