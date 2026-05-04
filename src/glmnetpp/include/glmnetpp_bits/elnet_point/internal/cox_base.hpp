#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/base.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/cox_adapter.hpp>

namespace glmnetpp {

/*
 * Base class for internal implementation of Cox elastic-net point solver.
 * This contains all the common interface and members across Cox versions.
 *
 * Cox regression uses IRLS (iteratively reweighted least squares):
 *   - Working response: z = eta + score / info
 *   - Working weights:  w = -diag_hessian (positive)
 *
 * State is stored in this base class following the binomial/poisson pattern:
 *   - r_: residuals (w * (z - eta))
 *   - w_: working weights
 *   - z_: working response
 *   - f_: linear predictor
 *
 * The CoxDevAdapter wraps coxdev library for IRLS integration while keeping
 * coxdev as a standalone, reusable library (usable by R and Python).
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalCoxBase
    : ElnetPointInternalNonLinearBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalNonLinearBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;

protected:
    using state_t = util::control_flow;
    using typename base_t::vec_t;
    using typename base_t::sp_mat_t;
    using mat_t = Eigen::Matrix<ValueType, Eigen::Dynamic, Eigen::Dynamic>;
    using ivec_t = Eigen::Matrix<IndexType, Eigen::Dynamic, 1>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    /*
     * Constructor for Cox point solver.
     *
     * @param thr       Convergence threshold
     * @param maxit     Maximum iterations
     * @param nx        Max number of variables to enter
     * @param nlp       Number of passes counter (reference)
     * @param ia        Active set indices (reference)
     * @param no        Number of observations
     * @param ni        Number of variables
     * @param dev0      Null deviance (reference, will be computed)
     * @param surv      Survival data structure (event times, status, etc.)
     * @param g         Linear predictor / offset (eta)
     * @param q         Original sample weights
     * @param vp        Penalty factors
     * @param cl        Box constraints
     * @param ju        Variable inclusion flags
     * @param int_param Internal parameters
     */
    template <class IAType
            , class SurvType
            , class GType
            , class QType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType
            >
    ElnetPointInternalCoxBase(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            index_t no,
            index_t ni,
            value_t& dev0,
            const SurvType& surv,
            const GType& g,
            const QType& q,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(thr, maxit, nx, nlp, false /* no intercept for Cox */, ia, dev0, vp, cl, ju)
        , a_(ni)
        , as_(ni)
        , xv_(ni)        // weighted variance of x columns
        , r_(no)         // residuals: w * (z - eta)
        , w_(no)         // working weights
        , z_(no)         // working response
        , f_(no)         // linear predictor
        , f_prev_(no)    // previous linear predictor (for IRLS convergence)
        , score_(no)     // score (gradient of log-likelihood)
        , score_prev_(no) // previous score (for IRLS convergence)
        , irls_thr_(thr) // IRLS convergence threshold (unscaled)
        , irls_change_(std::numeric_limits<value_t>::max())
        , q_(q.data(), q.size())
        , g_(g.data(), g.size())
        , cox_adapter_()
    {
        a_.setZero();
        as_.setZero();

        // Initialize adapter with survival data AND weights
        // The new coxdev requires weights at construction time
        // (used in preprocessing sort order for zero-weight handling)
        cox_adapter_.initialize(surv, q_);

        // Initialize linear predictor to offset
        f_ = g_;
    }

    // Cox has no intercept
    GLMNETPP_STRONG_INLINE auto intercept() const { return 0.0; }
    GLMNETPP_STRONG_INLINE auto beta(index_t k) const { return a_(k); }

    GLMNETPP_STRONG_INLINE
    void update_dlx(index_t k, value_t beta_diff) {
        base_t::update_dlx(beta_diff, xv_(k));
    }

    template <class PackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PackType& p)
    {
        base_t::compute_strong_map(
                this->abs_grad(), this->penalty(), this->strong_map(),
                p.beta, p.alm, p.alm0,
                [&](auto k) { return !this->is_excluded(k) || !this->exclusion()[k]; });
    }

    // Returns fraction of null deviance explained = 1 - current_dev / null_dev
    GLMNETPP_STRONG_INLINE value_t deviance() const {
        return 1.0 - cox_adapter_.deviance() / cox_adapter_.null_deviance();
    }

    GLMNETPP_STRONG_INLINE const auto& prediction() const { return f_; }

protected:
    GLMNETPP_STRONG_INLINE auto& beta(index_t k) { return a_(k); }
    GLMNETPP_STRONG_INLINE auto& resid() { return r_; }
    GLMNETPP_STRONG_INLINE const auto& resid() const { return r_; }
    GLMNETPP_STRONG_INLINE auto& x_var() { return xv_; }
    GLMNETPP_STRONG_INLINE const auto& x_var() const { return xv_; }
    GLMNETPP_STRONG_INLINE auto& weight() { return w_; }
    GLMNETPP_STRONG_INLINE const auto& weight() const { return w_; }
    GLMNETPP_STRONG_INLINE auto& working_response() { return z_; }
    GLMNETPP_STRONG_INLINE const auto& working_response() const { return z_; }
    GLMNETPP_STRONG_INLINE auto& linear_pred() { return f_; }
    GLMNETPP_STRONG_INLINE const auto& linear_pred() const { return f_; }
    GLMNETPP_STRONG_INLINE const auto& score() const { return score_; }
    GLMNETPP_STRONG_INLINE const auto& orig_weight() const { return q_; }
    GLMNETPP_STRONG_INLINE const auto& offset() const { return g_; }
    GLMNETPP_STRONG_INLINE auto& adapter() { return cox_adapter_; }
    GLMNETPP_STRONG_INLINE const auto& adapter() const { return cox_adapter_; }

    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t gk, value_t l1_regul, value_t l2_regul) {
        const auto& cl = this->endpts();
        base_t::update_beta(
                beta(k), gk, xv_(k), this->penalty()(k),
                cl(0,k), cl(1,k), l1_regul, l2_regul);
    }

    /*
     * Save old coefficients and update weighted variances.
     * Called at start of each WLS iteration.
     */
    template <class UpdateXVType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(UpdateXVType update_xv_f) {
        // Save old coefficients for convergence check
        std::for_each(this->active_begin(), this->active_end(),
                [&](auto k) { as_(k) = a_(k); });

        // Update weighted variance of each x column
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                update_xv_f,
                [&](auto k) { return this->is_excluded(k); });
    }

    /*
     * Recompute IRLS quantities using the adapter.
     * Called once per outer IRLS iteration.
     * Fills w_, z_, r_ from adapter computation.
     */
    /*
     * Save current state before IRLS recomputation.
     * Must be called BEFORE updating f_ in the derived class's update_irls().
     */
    GLMNETPP_STRONG_INLINE
    void save_irls_state() {
        f_prev_ = f_;
        score_prev_ = score_;
    }

    GLMNETPP_STRONG_INLINE
    void recompute_irls_quantities() {
        cox_adapter_.compute_irls_quantities(f_, q_, w_, z_, r_);
        cox_adapter_.get_score(score_);

        // Compute IRLS convergence measure (adelie-style):
        // |sum((score_new - score_old) * (eta_new - eta_old))|
        irls_change_ = std::abs((score_ - score_prev_).dot(f_ - f_prev_));
    }

    /*
     * Initialize the Cox point solver.
     * Computes null deviance and initial gradient.
     */
    template <class AbsGradFType>
    void construct(AbsGradFType abs_grad_f)
    {
        auto& dev0 = this->null_deviance();

        // Compute null deviance and initial IRLS quantities
        recompute_irls_quantities();
        cox_adapter_.set_null_deviance(cox_adapter_.deviance());
        dev0 = cox_adapter_.null_deviance();

        // Scale inner WLS threshold by null deviance (for coordinate descent)
        this->set_thresh(this->thresh() * dev0);

        // Initialize prev vectors for IRLS convergence check
        f_prev_ = f_;
        score_prev_ = score_;
        irls_change_ = std::numeric_limits<value_t>::max();

        // Compute initial absolute gradients for lambda_max
        gaussian_naive_t::update_abs_grad(this->abs_grad(), abs_grad_f,
                [&](index_t j) { return !this->exclusion()[j]; });
    }

    /*
     * Check IRLS convergence and update gradient.
     *
     * Convergence requires:
     * 1. Coefficient changes are small (checked by has_converged_irls)
     * 2. For unpenalized (lambda=0) cases: score (gradient) must also be small
     *
     * The additional score check for lambda=0 is critical for flat likelihood
     * surfaces (e.g., stratified Cox with ties) where the objective may plateau
     * while the gradient is non-zero. Without this check, the algorithm can
     * terminate prematurely before reaching the true optimum.
     */
    template <class GradFType>
    GLMNETPP_STRONG_INLINE
    state_t check_irls_convergence(
            value_t l1_regul,
            value_t l2_regul,
            GradFType grad_f)
    {
        // IRLS convergence: check if the quadratic approximation has stabilized
        // Uses adelie-style criterion: |sum((score_new - score_old) * (eta_new - eta_old))|
        // This measures change in the loss function and is naturally scaled.
        // Uses separate irls_thr_ (unscaled input thresh) to avoid the dev0-scaling
        // issue that caused 50-100x regression with coefficient-based criterion.
        bool ix = (irls_change_ <= irls_thr_);

        if (ix) {
            // For unpenalized (lambda ≈ 0) cases, also check score norm.
            // At the optimum, the score X_k^T * gradient should be ~0.
            // This is essential for convergence in flat likelihood regions
            // (e.g. lambda=0 with heavy ties).
            // Gate on BOTH l1 and l2 being zero: l1_regul = lambda*alpha is
            // also zero for pure ridge (alpha=0, lambda>0), but the ridge
            // optimum has X_k^T*score = lambda*beta_k != 0, so requiring the
            // score-norm to be small there makes the IRLS loop never converge.
            if (l1_regul < 1e-10 && l2_regul < 1e-10) {
                value_t max_active_score = 0.0;
                for (auto it = this->active_begin(); it != this->active_end(); ++it) {
                    auto k = *it;
                    max_active_score = std::max(max_active_score, std::abs(grad_f(k)));
                }

                // Score should be small relative to threshold.
                if (max_active_score > this->thresh()) {
                    // Coefficients appear converged but gradient is not small enough.
                    return state_t::noop_;
                }
            }

            // Check KKT for inactive variables
            auto skip_f = [&](auto k) { return !this->is_excluded(k) || !this->exclusion()[k]; };
            gaussian_naive_t::update_abs_grad(
                    this->abs_grad(), [&](index_t k) { return std::abs(grad_f(k)); }, skip_f);
            bool kkt_passed = gaussian_naive_t::check_kkt(
                    this->abs_grad(), this->penalty(), this->strong_map(),
                    l1_regul, skip_f);
            if (kkt_passed) return state_t::break_;
        }
        return state_t::noop_;
    }

private:
    // === State stored in base class (following binomial/poisson pattern) ===
    vec_t a_;                   // coefficients
    vec_t as_;                  // old coefficients (for convergence check)
    vec_t xv_;                  // weighted variance of columns of x
    vec_t r_;                   // residuals: w * (z - eta)
    vec_t w_;                   // working weights
    vec_t z_;                   // working response
    vec_t f_;                   // linear predictor
    vec_t f_prev_;              // previous linear predictor (for IRLS convergence)
    vec_t score_;               // score (gradient of log-likelihood) for KKT checks
    vec_t score_prev_;          // previous score (for IRLS convergence)
    value_t irls_thr_;          // IRLS convergence threshold (unscaled)
    value_t irls_change_;       // current IRLS convergence measure
    Eigen::Map<const vec_t> q_; // original sample weights
    Eigen::Map<const vec_t> g_; // offset / initial linear predictor

    // Adapter wrapping coxdev library
    CoxDevAdapter<value_t, index_t> cox_adapter_;
};

} // namespace glmnetpp
