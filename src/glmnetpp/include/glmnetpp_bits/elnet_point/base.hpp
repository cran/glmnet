#pragma once
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_point/traits.hpp>
#include <glmnetpp_bits/util/functional.hpp>
#include <glmnetpp_bits/elnet_point/decl.hpp>

namespace glmnetpp {

/*
 * Common CRTP base class for all point-solvers.
 * Only routines that __must__ be generated differently should be in here,
 * i.e. those that are dependent on the derived type.
 */
template <class ElnetPointDerived>
struct ElnetPointCRTPBase: 
    details::traits<ElnetPointDerived>::internal_t
{
private:
    using derived_t = ElnetPointDerived;
    using internal_t = typename details::traits<derived_t>::internal_t;

protected:
    using update_t = util::update_type;
    using typename internal_t::value_t;
    using typename internal_t::index_t;
    using state_t = util::control_flow;

    // Generate CRTP self()
    GLMNETPP_GENERATE_CRTP(derived_t)

    template <update_t upd, bool do_kkt, class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    std::pair<bool, bool> fit(const PointConfigPack& pack)
    {
        this->increment_passes(); 
        this->coord_desc_reset();
        util::if_else<upd == update_t::full>(
                [this, &pack]() {
                    this->for_each_with_skip(
                        this->all_begin(),
                        this->all_end(), 
                        [=, &pack](auto k) { this->self().template update<update_t::full>(k, pack); },
                        [=](auto k) { return this->is_excluded(k); }
                    );
                },
                [this, &pack]() {
                    this->for_each_with_skip(
                        this->active_begin(),
                        this->active_end(),
                        [=, &pack](auto k) { this->self().template update<update_t::partial>(k, pack); },
                        [](auto) { return false; } // no skip
                        );
                });
        this->update_intercept();
        if (this->has_converged()) { 
            return util::if_else<do_kkt>(
                    [this, &pack]() -> std::pair<bool, bool> { return {true, this->check_kkt(pack)}; },
                    []() -> std::pair<bool, bool> { return {true, true}; }
                    );
        }
        if (this->has_reached_max_passes()) { 
            throw util::maxit_reached_error();
        }
        return {false, false};
    }

    template <update_t upd, class PointPackType, class DiffType>
    GLMNETPP_STRONG_INLINE
    state_t update(index_t k, const PointPackType& pack, DiffType&& diff)
    {
        diff = this->beta(k);           // save old beta_k
        this->update_beta(k, pack);     // update new beta_k (assumes diff doesn't change)

        if (this->equal(diff, this->beta(k))) return state_t::continue_;

        // update active set stuff if full
        util::if_else<upd == update_t::full>(
                [=]() {
                    if (!this->is_active(k)) {
                        this->update_active(k);
                    }
                },
                []() {});

        diff = this->beta(k) - diff;    // new minus old beta_k
        this->update_dlx(k, diff);

        return state_t::noop_;
    }

public:
    using internal_t::internal_t;
};

/*
 * This is a CRTP base class for all non-linear GLM point-solvers.
 * All of these use IRLS-WLS implementation, so it is convenient capture the logic in one place.
 */
template <class ElnetPointDerived>
struct ElnetPointNonLinearCRTPBase: 
    ElnetPointCRTPBase<ElnetPointDerived>
{
private:
    using base_t = ElnetPointCRTPBase<ElnetPointDerived>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::update_t;
    using base_t::self;

public:
    using base_t::base_t;

    template <class PointConfigPack>
    void fit(const PointConfigPack& pack) { self().irls(pack); }

    /* 
     * Default version of coordinate-descent update logic that propagates update of coefficient difference.
     */
    template <update_t upd, class PointConfigPack, class DiffType>
    GLMNETPP_STRONG_INLINE
    state_t update(index_t k, const PointConfigPack& pack, DiffType&& diff)
    {
        state_t state = base_t::template update<upd>(k, pack, diff);
        if (state == state_t::continue_) return state_t::continue_;
        this->update_resid(k, diff);
        return state_t::noop_;
    }

    /* 
     * Default version of coordinate-descent update logic that should be called from the most derived object.
     */
    template <update_t upd, class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    void update(index_t k, const PointConfigPack& pack)
    {
        value_t diff = 0.0;
        self().template update<upd>(k, pack, diff);
    }

    /* 
     * Default version of WLS loop logic.
     */
    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    void wls(const PointConfigPack& pack) 
    {
        while (1) 
        {
            bool converged = false, _ = false;
            std::tie(converged, _) = base_t::template fit<update_t::full, false>(pack);
            if (converged) break;

            // partial fit
            while (1) {
                bool converged = false, _ = false;
                std::tie(converged, _) = base_t::template fit<update_t::partial, false>(pack);
                if (converged) break; 
            }
        }
    }

    /* 
     * Default version of IRLS loop logic.
     */
    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    void irls(const PointConfigPack& pack) 
    {
        this->initialize(pack);

        while (1) {
            if (this->has_reached_max_passes()) {
                throw util::maxit_reached_error();
            }
            this->setup_wls(pack);
            base_t::self().wls(pack);
            state_t state = this->update_irls(pack);
            if (state == state_t::break_) break;
        }
    }
};

} // namespace glmnetpp
