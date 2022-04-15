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
    util::mode_type<util::glm_type::binomial>::multi_class,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalBinomialMultiClassBase<ValueType, IndexType, BoolType>
        , SpElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialMultiClassBase<ValueType, IndexType, BoolType>;
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
            , class ISType
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
            GType& g,
            value_t& dev0,
            const XType& X,
            const YType& y,
            const WType& w,
            const XBType& xb,
            const XSType& xs,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            ISType& is,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, g, dev0, y, w, vp, cl, ju, is, int_param)
        , sp_base_t(X, xb, xs)
    {
        base_t::construct(
                [&](index_t j) { return sp_base_t::compute_xv(j, this->weight()); },
                [&](index_t ic) { initialize_resid(ic); },
                [&](index_t j) { return sp_base_t::compute_grad(j, this->resid(), this->new_weight()); });
    }

    using base_t::check_kkt;
    using base_t::update_dlx;
    using base_t::for_each_with_skip;

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
    void update_beta(index_t k, const PointPackType& pack) {
        base_t::update_beta(k, sp_base_t::compute_grad(k, this->resid(), this->new_weight()), pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    state_t setup_wls(index_t ic) {
        state_t state = base_t::setup_wls(ic);
        if (state != state_t::noop_) return state;
        const auto& v = this->new_weight();
        auto xmz = this->new_weight_sum();
        auto& xv_ic = this->curr_x_var();
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                [&](index_t j) { 
                    sp_base_t::update_with_new_weights(j, v, this->optimization_type(), xmz, xv_ic[j]);
                },
                [&](index_t j) { return this->is_excluded(j); });

        sp_base_t::update_shifts(this->resid().sum());
        return state_t::noop_;
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        value_t b0 = 0.0;
        return base_t::update_irls(
                pack.elastic_prop(), pack.l1_regul(),
                [this, &b0](index_t l, value_t s, auto& buff) {
                    this->sp_base_t::update_prediction(l, s, buff, b0);
                },
                [&b0](auto& buff) { buff.array() = (buff.array() + b0).exp(); },
                [&](index_t ic) { initialize_resid(ic); },
                [this](index_t k) { return sp_base_t::compute_grad(k, this->resid(), this->new_weight()); });
    }

    GLMNETPP_STRONG_INLINE
    void update_irls_class()
    {
        base_t::has_converged_irls_class();
        base_t::update_irls_class(
                [&](auto& buff) {
                    value_t b0 = 0.0;
                    std::for_each(this->active_begin(), this->active_end(),
                        [&](index_t k) { sp_base_t::update_prediction(k, -this->beta(k), buff, b0); });
                    buff.array() += b0;
                });
    }

private:
    GLMNETPP_STRONG_INLINE
    void initialize_resid(index_t ic) {
        auto& y = this->y();
        auto& v = this->new_weight();
        auto& r = this->resid();
        auto& q = this->q();
        auto& sxp = this->sxp();
        v.array() = q.col(ic).array()/sxp.array(); 
        base_t::initialize_resid(r, y.col(ic), v);
        v.array() = this->weight().array() * v.array() * (1.0 - v.array());
        sp_base_t::update_shifts(r.sum());
    }
};

} // namespace glmnetpp

