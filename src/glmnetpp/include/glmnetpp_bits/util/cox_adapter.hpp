#pragma once
#include <cmath>
#include <memory>
#include <Eigen/Core>
#include <coxdev.h>

namespace glmnetpp {

/**
 * Lightweight survival data structure for passing through glmnetpp layers.
 * Constructed at the R/C++ boundary and passed to driver -> path -> point.
 */
template <class ValueType, class IndexType>
struct CoxSurvivalData {
    using value_t = ValueType;
    using index_t = IndexType;
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using ivec_t = Eigen::Matrix<index_t, Eigen::Dynamic, 1>;

    vec_t start;   // start times (for left-truncation)
    vec_t stop;    // stop/event times
    ivec_t status; // event indicator (1=event, 0=censored)
    ivec_t strata; // strata labels (empty = single stratum)
    bool efron;    // true for Efron, false for Breslow

    template <class StartType, class StopType, class StatusType>
    CoxSurvivalData(const StartType& start_, const StopType& stop_,
                    const StatusType& status_, bool efron_)
        : start(start_), stop(stop_), status(status_), strata(), efron(efron_)
    {}

    template <class StartType, class StopType, class StatusType, class StrataType>
    CoxSurvivalData(const StartType& start_, const StopType& stop_,
                    const StatusType& status_, const StrataType& strata_, bool efron_)
        : start(start_), stop(stop_), status(status_), strata(strata_), efron(efron_)
    {}
};

/**
 * Adapter that wraps the coxdev library for glmnetpp IRLS integration.
 *
 * This adapter:
 * - Owns a StratifiedCoxDeviance (or CoxDeviance) from the coxdev library
 * - Provides compute_irls_quantities() that fills OUTPUT buffers owned by caller
 * - Keeps coxdev as a standalone, reusable library
 *
 * The expensive preprocessing happens once at initialization (constructor).
 * Each IRLS iteration reuses the internal buffers managed by coxdev.
 */
template <class ValueType, class IndexType>
struct CoxDevAdapter {
    using value_t = ValueType;
    using index_t = IndexType;
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using ivec_t = Eigen::Matrix<index_t, Eigen::Dynamic, 1>;

private:
    // Owns either a CoxDeviance or StratifiedCoxDeviance
    std::unique_ptr<CoxDeviance> cox_dev_;
    std::unique_ptr<StratifiedCoxDeviance> strat_cox_dev_;
    bool is_stratified_ = false;

    // Pre-allocated intermediate buffers
    vec_t score_buffer_;  // score = d(loglik)/d(eta)

    int n_total_ = 0;
    value_t deviance_ = 0.0;
    value_t null_deviance_ = 0.0;
    bool initialized_ = false;

public:
    CoxDevAdapter() = default;

    /**
     * Initialize from survival data and weights.
     * Performs preprocessing (O(n log n)) and allocates internal buffers.
     * Call once at construction time.
     *
     * Note: The new coxdev requires weights at construction time
     * (used in preprocessing sort order for zero-weight handling).
     */
    template <class SurvType, class WeightType>
    void initialize(const SurvType& surv, const WeightType& weights) {
        n_total_ = surv.status.size();
        score_buffer_.resize(n_total_);

        // Convert to Eigen types expected by coxdev
        Eigen::VectorXd start_d = surv.start.template cast<double>();
        Eigen::VectorXd stop_d = surv.stop.template cast<double>();
        Eigen::VectorXi status_i = surv.status.template cast<int>();
        Eigen::VectorXd weight_d = weights.template cast<double>();

        if (surv.strata.size() > 0) {
            Eigen::VectorXi strata_i = surv.strata.template cast<int>();
            strat_cox_dev_ = std::make_unique<StratifiedCoxDeviance>(
                start_d, stop_d, status_i, strata_i, weight_d, surv.efron);
            is_stratified_ = true;
        } else {
            cox_dev_ = std::make_unique<CoxDeviance>(
                start_d, stop_d, status_i, weight_d, surv.efron);
            is_stratified_ = false;
        }

        initialized_ = true;
    }

    /**
     * Compute IRLS quantities and write to OUTPUT buffers.
     * Called once per outer IRLS iteration.
     *
     * @param eta                    Linear predictor (input)
     * @param weights                Sample weights (input)
     * @param working_weights_out    OUTPUT: working weights for WLS (w = -diag_hessian / 2)
     * @param working_response_out   OUTPUT: working response (z = eta + score/info)
     * @param residuals_out          OUTPUT: residuals (r = w * (z - eta))
     * @return Current deviance
     */
    template <class EtaType, class WeightType>
    value_t compute_irls_quantities(
            const EtaType& eta,
            const WeightType& weights,
            vec_t& working_weights_out,
            vec_t& working_response_out,
            vec_t& residuals_out)
    {
        // Convert to Eigen::VectorXd for coxdev
        Eigen::VectorXd eta_d = eta.template cast<double>();
        Eigen::VectorXd w_d = weights.template cast<double>();

        // Call coxdev
        if (is_stratified_) {
            deviance_ = strat_cox_dev_->compute_deviance(eta_d, w_d);
        } else {
            deviance_ = cox_dev_->compute_deviance(eta_d, w_d);
        }

        // Get gradient and diagonal Hessian from coxdev
        // coxdev returns: gradient = -2 * d(loglik)/d(eta)
        //                 diag_hessian = -2 * d²(loglik)/d(eta)²
        const Eigen::VectorXd& grad = is_stratified_
            ? strat_cox_dev_->get_gradient()
            : cox_dev_->get_gradient();
        const Eigen::VectorXd& diag_h = is_stratified_
            ? strat_cox_dev_->get_diag_hessian()
            : cox_dev_->get_diag_hessian();

        // Convert to glmnetpp working quantities:
        //   score = d(loglik)/d(eta) = -gradient / 2
        //   working_weights = -diag_hessian / 2 (positive for minimization)
        int n = eta.size();
        for (int i = 0; i < n; ++i) {
            score_buffer_(i) = static_cast<value_t>(-grad(i) * 0.5);
            working_weights_out(i) = static_cast<value_t>(diag_h(i) * 0.5);
        }

        // Working response: z = eta + score / info
        for (int i = 0; i < n; ++i) {
            if (std::abs(working_weights_out(i)) > 1e-10) {
                working_response_out(i) = eta(i) + score_buffer_(i) / working_weights_out(i);
            } else {
                working_response_out(i) = eta(i);
            }
        }

        // Residuals: r = w * (z - eta)
        residuals_out.array() = working_weights_out.array() * (working_response_out - eta).array();

        return static_cast<value_t>(deviance_);
    }

    /**
     * Compute null deviance (deviance at eta=0 or eta=offset).
     * Should be called once after initialization to set baseline.
     */
    template <class EtaType, class WeightType>
    value_t compute_null_deviance(const EtaType& eta, const WeightType& weights) {
        vec_t dummy_w(eta.size()), dummy_z(eta.size()), dummy_r(eta.size());
        null_deviance_ = compute_irls_quantities(eta, weights, dummy_w, dummy_z, dummy_r);
        return null_deviance_;
    }

    /**
     * Get the score (gradient of log-likelihood) from the last computation.
     * Useful for lambda_max calculation and KKT checking.
     * Note: After compute_irls_quantities, score_buffer_ already contains the score.
     */
    void get_score(vec_t& score_out) const {
        score_out = score_buffer_;
    }

    // Accessors
    value_t deviance() const { return static_cast<value_t>(deviance_); }
    value_t null_deviance() const { return static_cast<value_t>(null_deviance_); }
    void set_null_deviance(value_t d) { null_deviance_ = d; }

    int n_total() const { return n_total_; }
    bool is_stratified() const { return is_stratified_; }
    bool is_initialized() const { return initialized_; }
};

} // namespace glmnetpp
