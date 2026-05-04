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
#include <glmnetpp_bits/elnet_point/internal/cox_base.hpp>

namespace glmnetpp {

/*
 * Dense elastic-net point solver for Cox naive method.
 *
 * Uses CoxDevAdapter in base class to wrap coxdev library.
 * State is stored in base class following the binomial/poisson pattern:
 *   - r_: residuals
 *   - w_: working weights
 *   - z_: working response
 *   - f_: linear predictor
 *
 * Single-stratum models are treated as the special case of one stratum.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
    util::glm_type::cox,
    util::mode_type<util::glm_type::cox>::naive,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalCoxBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalCoxBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
        ValueType, IndexType, BoolType>;
    using typename base_t::state_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class XType
            , class SurvType
            , class GType
            , class QType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternal(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            value_t& dev0,
            const XType& X,
            const SurvType& surv,
            const GType& g,
            const QType& q,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(thr, maxit, nx, nlp, ia, X.rows(), X.cols(), dev0, surv, g, q, vp, cl, ju, int_param)
        , X_(X.data(), X.rows(), X.cols())
    {
        // Base class constructor:
        // - Initializes adapter with survival data (preprocessing + workspaces)
        // - Initializes linear predictor to offset

        // Compute null deviance and initial gradient using base class construct()
        base_t::construct([&](index_t j) { return compute_abs_grad(j); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        auto gk = compute_grad(k);
        base_t::update_beta(k, gk, pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        // Update residuals incrementally: r -= delta * w * x_k
        this->resid().array() -= beta_diff * this->weight().array() * X_.col(k).array();
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        // Cox regression has no intercept
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(const PointPackType&) {
        // Compute weighted variance of each x column: v(k) = X_k^T * diag(w) * X_k
        const auto& w = this->weight();
        base_t::setup_wls([&](auto k) {
            this->x_var()(k) = compute_xv(X_.col(k), w);
        });
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack)
    {
        // Save previous state BEFORE updating f
        this->save_irls_state();

        // Update linear predictor: f = offset + X * beta
        auto& f = this->linear_pred();
        f = this->offset();
        std::for_each(this->active_begin(), this->active_end(),
                [&](index_t k) { f += this->beta(k) * X_.col(k); });

        // Recompute IRLS quantities at new linear predictor
        this->recompute_irls_quantities();

        // Check convergence using loss-change criterion (adelie-style)
        value_t l1_regul = pack.l1_regul();
        value_t l2_regul = pack.l2_regul();
        auto grad_f = [&](index_t k) { return compute_grad_from_score(k); };
        return base_t::check_irls_convergence(l1_regul, l2_regul, grad_f);
    }

private:
    using typename base_t::mat_t;
    using typename base_t::vec_t;
    using typename base_t::ivec_t;

    template <class XColType, class WType>
    GLMNETPP_STRONG_INLINE
    static value_t
    compute_xv(const XColType& x, const WType& w) {
        return w.dot(x.array().square().matrix());
    }

    // Gradient for coordinate descent update (uses current WLS residual)
    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t j) const {
        return X_.col(j).dot(this->resid());
    }

    // Gradient from Cox score function (used for lambda_max and KKT checking)
    GLMNETPP_STRONG_INLINE
    value_t compute_grad_from_score(index_t j) const {
        return X_.col(j).dot(this->score());
    }

    // Absolute gradient for lambda_max computation: |X_j^T * score|
    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t j) const {
        return std::abs(compute_grad_from_score(j));
    }

    Eigen::Map<const mat_t> X_;
};

} // namespace glmnetpp
