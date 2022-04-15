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
 * Dense elastic-net point solver for Binomial multi-class method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
    util::glm_type::binomial,
    util::mode_type<util::glm_type::binomial>::multi_class,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalBinomialMultiClassBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialMultiClassBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
        ValueType, IndexType, BoolType>;
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
            , class VPType
            , class CLType
            , class JUType
            , class ISType
            , class IntParamType>
    ElnetPointInternal(
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
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            ISType& is,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, g, dev0, y, w, vp, cl, ju, is, int_param)
        , X_(X.data(), X.rows(), X.cols())
    {
        base_t::construct(
                [&](index_t j) { return compute_xv(j, this->weight()); },
                [&](index_t ic) { base_t::initialize_resid(ic); },
                [&](index_t j) { return compute_grad(j); });
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() { base_t::update_intercept(this->resid().sum()); }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        gaussian_naive_t::update_resid(
                this->resid(), beta_diff, 
                (this->new_weight().array() * X_.col(k).array()).matrix());
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        base_t::update_beta(k, compute_grad(k), pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    state_t setup_wls(index_t ic) {
        state_t state = base_t::setup_wls(ic);
        if (state != state_t::noop_) return state;
        if (!this->optimization_type()) {
            auto& xv_ic = this->curr_x_var();
            base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                    [&](index_t j) { xv_ic(j) = this->compute_xv(j, this->new_weight()); },
                    [&](index_t j) { return this->is_excluded(j); });
        }
        return state_t::noop_;
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        return base_t::update_irls(pack.elastic_prop(), pack.l1_regul(), 
                [&](index_t l, value_t s, auto& buff) { buff -= s * X_.col(l); },
                [&](auto& buff) { buff.array() = buff.array().exp(); },
                [&](index_t ic) { base_t::initialize_resid(ic); }, 
                [&](index_t k) { return this->compute_grad(k); });
    }

    GLMNETPP_STRONG_INLINE
    void update_irls_class()
    {
        base_t::has_converged_irls_class();
        base_t::update_irls_class(
                [&](auto& buff) {
                    std::for_each(this->active_begin(), this->active_end(),
                        [&](index_t k) { buff += this->beta(k) * X_.col(k); });
                });
    }

private:
    using typename base_t::mat_t;

    template <class WType>
    GLMNETPP_STRONG_INLINE
    value_t compute_xv(index_t j, const WType& w) const { 
        return w.dot(X_.col(j).array().square().matrix());
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t k) const {
        return base_t::compute_grad(this->resid(), X_.col(k));
    }

    Eigen::Map<const mat_t> X_;
};

} // namespace glmnetpp
