#pragma once
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
    using typename base_t::update_type;
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using base_t::partial_fit;
    using base_t::update;

public:
    using base_t::base_t;

    void partial_fit(index_t m, value_t ab, value_t dem)
    {
        this->compress_active();

        base_t::partial_fit(m, ab, dem);

        // update gradient to leave invariant
        this->update_compressed_active();
        this->update_grad_compressed_active();
    }

    template <update_type upd>
    void update(index_t k, value_t ab, value_t dem)
    {
        value_t beta_diff = 0;

        auto state = base_t::template update<upd>(k, ab, dem, beta_diff);
        if (state == state_t::continue_) return;

        // update gradient
        util::if_else<upd == update_type::full>(
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
