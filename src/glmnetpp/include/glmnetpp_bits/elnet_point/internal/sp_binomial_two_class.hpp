#pragma once
#include <cmath>
#include <utility>
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/iterator/counting_iterator.hpp>
#include <glmnetpp_bits/util/iterator/one_to_zero_iterator.hpp>
#include <glmnetpp_bits/elnet_point/internal/binomial_base.hpp>

namespace glmnetpp {

/*
 * Sparse elastic-net point solver for Binomial multi-class method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternal<
    util::glm_type::binomial,
    util::mode_type<util::glm_type::binomial>::two_class,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalBinomialTwoClassBase<ValueType, IndexType, BoolType>
        , SpElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialTwoClassBase<ValueType, IndexType, BoolType>;
    using sp_base_t = SpElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>;
    using typename base_t::state_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class GType
            , class XType
            , class YType
            , class WType
            , class XBType
            , class XSType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    SpElnetPointInternal(
            bool isd,
            bool intr,
            index_t kopt,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const GType& g,
            value_t& dev0,
            const XType& X,
            const YType& y,
            const WType& w,
            const XBType& xb,
            const XSType& xs,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, g, dev0, y, w, vp, cl, ju, int_param)
        , sp_base_t(X, xb, xs)
        , sc_(X.rows())
    {
        this->construct(
                [&](index_t j) { return sp_base_t::compute_xv(j, this->weight()); },
                [&](index_t j) { return sp_base_t::compute_grad(j, this->resid(), this->new_weight()); });
        sp_base_t::update_shifts(this->resid().sum());
    }

    using base_t::check_kkt;
    using base_t::update_dlx;
    using base_t::for_each_with_skip;

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        auto gk = sp_base_t::compute_grad(k, this->resid(), this->new_weight());
        base_t::update_beta(k, gk, pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    void update_active(index_t k) {
        base_t::update_active(k);
        sp_base_t::update_active(k, this->new_weight());
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        auto d = base_t::update_intercept(this->sum_weighted_resid());
        sp_base_t::update_intercept(d, this->new_weight_sum());
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        sp_base_t::update_resid(
                k, this->resid(), beta_diff, this->new_weight(), this->new_weight_sum());
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(const PointPackType&) {
        const auto& ixx = this->strong_map();
        const auto& v = this->new_weight();
        auto& xv = this->x_var();
        const auto& xmz = this->new_weight_sum();

        base_t::setup_wls();
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                [&](index_t j) { 
                    sp_base_t::update_with_new_weights(j, v, this->optimization_type(), xmz, xv(j));
                },
                [&](index_t j) { return !ixx[j]; });
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        sc_.array() = this->intercept(); 
        auto b0=0.0;
        std::for_each(this->active_begin(), this->active_end(),
                [&](index_t l) {
                    sp_base_t::update_prediction(l, -this->beta(l), sc_, b0);
                });
        sc_.array() += b0;

        auto predict_f = [&](index_t i) { return sc_(i) + this->offset()(i); };
        state_t state = base_t::update_irls_invariants(predict_f);
        if (state == state_t::break_) return state_t::break_;

        // update further invariants
        sp_base_t::update_shifts(this->resid().sum());

        auto grad_f = [&](index_t k) { return sp_base_t::compute_grad(k, this->resid(), this->new_weight()); };
        return base_t::update_irls_strong_set(grad_f, pack.l1_regul());
    }

private:
    using typename base_t::vec_t;
    
    vec_t sc_;  // buffer for temporary storage for IRLS
};

} // namespace glmnetpp
