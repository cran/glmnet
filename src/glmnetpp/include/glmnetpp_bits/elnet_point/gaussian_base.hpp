#pragma once
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_point/base.hpp>

namespace glmnetpp {

template <class ElnetPointGaussianDerived>
struct ElnetPointGaussianBase : 
    ElnetPointCRTPBase<ElnetPointGaussianDerived>
{
private:
    using base_t = ElnetPointCRTPBase<ElnetPointGaussianDerived>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::update_t;
    using typename base_t::state_t;
    using base_t::self;

public:
    using base_t::base_t;

    template <class PointConfigPack>
    void fit(const PointConfigPack& pack)
    {
        this->initialize(pack);

        if (this->is_warm_ever()) {
            self().partial_fit(pack);
        }

        while (1) {
            bool converged_kkt = this->initial_fit(
                    [&]() { return base_t::template fit<update_t::full, true>(pack); }
                    );
            if (converged_kkt) return;
            self().partial_fit(pack);
        }
    }

protected:

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void partial_fit(const PointPackType& pack)
    {
        this->set_warm_ever();

        // fit on partial subset
        while (1) {
            bool converged = false, _ = false;
            std::tie(converged, _) = base_t::template fit<update_t::partial, false>(pack);
            if (converged) break;
        }
    }

    template <update_t upd, class PointPackType, class DiffType>
    GLMNETPP_STRONG_INLINE
    state_t update(index_t k, const PointPackType& pack, DiffType&& diff)
    {
        state_t state = base_t::template update<upd>(k, pack, diff);
        if (state == state_t::continue_) return state_t::continue_;
        this->update_rsq(k, diff);
        return state_t::noop_;
    }
};

} // namespace glmnetpp
