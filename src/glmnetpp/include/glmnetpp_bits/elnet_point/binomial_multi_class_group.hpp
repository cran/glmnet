#pragma once
#include <glmnetpp_bits/elnet_point/decl.hpp>
#include <glmnetpp_bits/elnet_point/binomial_base.hpp>

namespace glmnetpp {

template <class ElnetPointInternalPolicy>
struct ElnetPoint<
    util::glm_type::binomial, 
    util::mode_type<util::glm_type::binomial>::multi_class_group,
    ElnetPointInternalPolicy>
        : ElnetPointBinomialBase<
            ElnetPoint<
                util::glm_type::binomial, 
                util::mode_type<util::glm_type::binomial>::multi_class_group,
                ElnetPointInternalPolicy> >
{
private:
    using base_t = ElnetPointBinomialBase<
        ElnetPoint<util::glm_type::binomial, 
                   util::mode_type<util::glm_type::binomial>::multi_class_group,
                   ElnetPointInternalPolicy> >;
    using typename base_t::index_t;
    using typename base_t::update_t;

public:
    using base_t::base_t;

    template <update_t upd, class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    void update(index_t k, const PointConfigPack& pack)
    {
        auto&& diff = this->beta_buffer();
        base_t::template update<upd>(k, pack, diff);
    }
};

} // namespace glmnetpp
