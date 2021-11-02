#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/traits.hpp>
#include <glmnetpp_bits/util/functional.hpp>

namespace glmnetpp {

template <class ElnetPointGaussianDerived>
struct ElnetPointGaussianBase : 
    protected details::traits<ElnetPointGaussianDerived>::internal_t
{
private:
    using derived_t = ElnetPointGaussianDerived;
    using internal_t = typename details::traits<derived_t>::internal_t;

    derived_t& self() { return static_cast<derived_t&>(*this); }
    const derived_t& self() const { return static_cast<const derived_t&>(*this); }

protected:
    using typename internal_t::value_t;
    using typename internal_t::index_t;

public:
    using internal_t::internal_t;
    using internal_t::beta;
    using internal_t::rsq;
    using internal_t::n_active;

    template <class PointConfigPack>
    void fit(const PointConfigPack& pack)
    {
        auto m = pack.m;
        auto ab = pack.ab;
        auto dem = pack.dem;
        auto alm0 = pack.alm0;
        auto alm = pack.alm;
        auto beta = pack.beta;

        this->reset_warm();

        this->initialize(beta, alm, alm0);

        if (this->is_warm_ever()) {
            self().partial_fit(m, ab, dem);
        }

        while (1) {

            auto initial_fit_internal = [=]() -> std::pair<bool, bool>{
                ++this->passes();
                this->coord_desc_reset();
                this->coord_desc(
                    this->all_begin(),
                    this->all_end(), 
                    [=](auto k) { self().template update<update_type::full>(k, ab, dem); },
                    [=](auto k) { return this->is_excluded(k); }
                );

                if (this->has_converged()) { 
                    return {true, this->check_kkt(ab)};
                }
                if (this->has_reached_max_passes()) { 
                    throw util::maxit_reached_error();
                }
                return {false, false};
            };

            bool converged_kkt = this->initial_fit(initial_fit_internal);

            if (converged_kkt) return;

            self().partial_fit(m, ab, dem);
        }
    }

protected:
    using state_t = util::control_flow;

    void partial_fit(index_t m, value_t ab, value_t dem)
    {
        this->set_warm_ever();

        // fit on partial subset
        while (1) {
            ++this->passes();
            this->coord_desc_reset();
            this->coord_desc(
                    this->active_begin(),
                    this->active_end(),
                    [=](auto k) { self().template update<update_type::partial>(k, ab, dem); },
                    [](auto) { return false; } // no skip
                    );
            if (this->has_converged()) break; 
            if (this->has_reached_max_passes()) {
                throw util::maxit_reached_error();
            }
        }

        this->set_warm();
    }

    enum class update_type
    {
        full,
        partial
    };

    template <update_type upd>
    state_t update(index_t k, value_t ab, value_t dem, value_t& diff)
    {
        auto ak = this->beta(k);
         
        this->update_beta(k, ab, dem);
        
        if (this->beta(k) == ak) return state_t::continue_;

        // update active set stuff if full
        util::if_else<upd == update_type::full>(
                [=]() {
                    if (!this->is_active(k)) {
                        this->update_active(k);
                    }
                },
                []() {});


        diff = this->beta(k) - ak; 
        this->update_rsq(k, diff);
        this->update_dlx(k, diff);

        return state_t::noop_;
    }
};

} // namespace glmnetpp
