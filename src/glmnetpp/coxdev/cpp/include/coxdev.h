#include <cstddef>
#include <limits>
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <cstdint>
#include <Eigen/Dense>

#ifdef PY_INTERFACE
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#ifdef R_INTERFACE
#include <RcppEigen.h>
#endif

#define EIGEN_REF Eigen::Ref
#define MAKE_MAP_Xd(y) Eigen::Map<Eigen::VectorXd>((y).data(), (y).size())
#define MAKE_MAP_Xi(y) Eigen::Map<Eigen::VectorXi>((y).data(), (y).size())

#if defined(PY_INTERFACE)
#define ERROR_MSG(x) throw std::runtime_error(x)
#elif defined(R_INTERFACE)
#define ERROR_MSG(x) Rcpp::stop(x)
#else
#define ERROR_MSG(x) throw std::runtime_error(x)
#endif

// The Breslow (or, with Efron's ties, Efron's) estimate of the baseline
// hazard at the linear predictor of the last compute_deviance: one entry per
// distinct failure time with positive weight, in increasing time within each
// stratum. hazard is the increment at that time and cumulative_hazard its
// running sum within the stratum, both for the linear predictor as given
// (not centered): the cumulative hazard of subject i at time t is
// cumulative_hazard(t) * exp(linear_predictor(i)). stratum is the stratum
// label (0 without strata).
struct BaselineHazard {
    Eigen::VectorXd time;
    Eigen::VectorXd hazard;
    Eigen::VectorXd cumulative_hazard;
    Eigen::VectorXi stratum;
};

// Harrell's concordance index (C-index) as a ratio: numerator / denominator,
// NaN when no pair is comparable (denominator 0). See concordance below.
struct Concordance {
    double numerator;
    double denominator;
    double value() const {
        return denominator > 0 ? numerator / denominator
                               : std::numeric_limits<double>::quiet_NaN();
    }
};

struct CoxContext {
    Eigen::Ref<const Eigen::VectorXi> event_order;
    Eigen::Ref<const Eigen::VectorXi> start_order;
    Eigen::Ref<const Eigen::VectorXi> status;
    Eigen::Ref<const Eigen::VectorXi> first;
    Eigen::Ref<const Eigen::VectorXi> last;
    Eigen::Ref<const Eigen::VectorXd> scaling;
    Eigen::Ref<const Eigen::VectorXi> event_map;
    Eigen::Ref<const Eigen::VectorXi> start_map;
    bool have_start_times;
    bool efron;
};

class CoxDeviance {
public:
    CoxDeviance(const Eigen::VectorXd& start,
                const Eigen::VectorXd& event,
                const Eigen::VectorXi& status,
                const Eigen::VectorXd& weight,
                bool efron);

    double compute_deviance(const Eigen::VectorXd& eta,
                            const Eigen::VectorXd& sample_weight);

    void compute_hessian_matvec(const Eigen::VectorXd& arg,
                                Eigen::VectorXd& out);

    // Draws of the score from the model at the linear predictor and weights of
    // the last compute_deviance (or, for "wild", from its score residuals). Each column r of the (n x n_draws) result is in native
    // order and on the scale of the score -gradient / 2, so that X^T r is a draw
    // of the score of the coefficients of a design X. The model-based methods
    // draw failures from the model and return their (Breslow) score: at each
    // event time, the sum over the failures of x_j minus its expectation.
    // method is one of
    //   "multinomial"  the d_k failures at each event time drawn with replacement
    //                  from the risk set, with probabilities proportional to
    //                  w exp(eta); the covariance is the Breslow information
    //   "bernoulli"    each subject at risk fails independently with
    //                  probability min(1, d_k p_k)
    //   "fisher"       the set of d_k failures drawn from Fisher's noncentral
    //                  hypergeometric distribution (independent Bernoullis given
    //                  d_k failures), centered at its exact mean: the score of
    //                  the exact discrete partial likelihood (Cox, 1972), by the
    //                  recursion of Gail, Lubin and Rubinstein (1981, Biometrika
    //                  68, 703-707); see the README for full references
    //   "wild"         sum_i G_i U_i over the score residuals U_i, with
    //                  independent standard normal G_i
    // The model-based methods (all but "wild") need equal positive weights.
    // Draws depend only on seed, not on the platform.
    Eigen::MatrixXd sample_score(const std::string& method,
                                   int n_draws,
                                   std::uint64_t seed);

    // The wild draw for given multipliers G (native order): r with
    // X^T r = sum_i G_i U_i, U_i the score residual of subject i (centered in
    // each risk set, with Efron's weights for ties if efron).
    Eigen::VectorXd wild_score(const Eigen::VectorXd& multipliers);

    // The score residuals of one covariate x (n-vector, native order) at the
    // linear predictor and weights of the last compute_deviance: entry i is
    // the contribution of subject i to the score x^T (-gradient / 2), centered
    // within each risk set (with Efron's weights for ties if efron). They are
    // linear in x, u = L x, and wild_score is the adjoint: G^T L x =
    // x^T wild_score(G), so they sum to the score. Applied to each column of a
    // design X they give the n x p score residuals U, whose crossproduct is
    // the middle of the robust (sandwich) variance.
    Eigen::VectorXd score_residuals(const Eigen::VectorXd& x) const;

    // The baseline hazard at the last compute_deviance; see BaselineHazard.
    BaselineHazard baseline_hazard() const;

    // Harrell's concordance index of a predictor pred (n-vector, native
    // order; a higher pred means a higher risk) for these data, with
    // weights weight; see the free function concordance. Does not need
    // compute_deviance.
    Concordance concordance(const Eigen::VectorXd& pred,
                            const Eigen::VectorXd& weight) const;

    // The Schoenfeld residuals of one covariate x (n-vector, native order) at
    // the last compute_deviance: for a failure i with positive weight, x_i
    // minus the mean of x over its risk set, with weights w exp(eta); with
    // Efron's ties, the mean is the average over the tied failures of their
    // Efron means, as in survival (Schoenfeld, 1982, Biometrika 69, 239-241;
    // Grambsch and Therneau, 1994, Biometrika 81, 515-526). 0 for the other
    // observations. Weighted by
    // w, they sum to the score x^T (-gradient / 2). Applied to each column
    // of a design X they give the Schoenfeld residuals of X.
    Eigen::VectorXd schoenfeld_residuals(const Eigen::VectorXd& x) const;

    // The risk-set variance of one covariate z (n-vector, native order) at
    // the last compute_deviance: for a failure i with positive weight, the
    // variance of z over its risk set with weights w exp(eta), the V_z of the
    // scaled Schoenfeld residuals (Grambsch and Therneau, 1994, Biometrika
    // 81, 515-526); with Efron's ties, the average over the tied failures of
    // their Efron variances, as for schoenfeld_residuals. 0 for the other
    // observations. Weighted by w, they sum to the information z^T I z.
    // Covariances of two covariates follow by polarization:
    // (V(y + z) - V(y - z)) / 4.
    Eigen::VectorXd risk_set_variance(const Eigen::VectorXd& z) const;

    // The expected number of events of each observation (native order) at
    // the last compute_deviance: exp(eta_i) times the increase of the
    // baseline cumulative hazard over (start_i, event_i], where with Efron's
    // ties a failure counts at its own failure time with Efron's weights.
    // As survival's predict(type = "expected"): status - expected events is
    // the martingale residual (Therneau, Grambsch and Fleming, 1990,
    // Biometrika 77, 147-160), and exp(-expected events) the probability of
    // surviving to event_i, given survival to start_i.
    Eigen::VectorXd expected_events() const;

    // Getters for buffers (for result construction or debugging)
    const Eigen::VectorXd& get_gradient() const { return grad_buffer; }
    const Eigen::VectorXd& get_diag_hessian() const { return diag_hessian_buffer; }
    const Eigen::VectorXd& get_linear_predictor() const { return linear_predictor; }
    const Eigen::VectorXd& get_sample_weight() const { return sample_weight; }
    double get_loglik_sat() const { return loglik_sat; }

    // Getters for preprocessed data
    const Eigen::VectorXi& get_event_order() const { return event_order; }
    const Eigen::VectorXi& get_start_order() const { return start_order; }
    const Eigen::VectorXi& get_first() const { return _first; }
    const Eigen::VectorXi& get_last() const { return _last; }
    const Eigen::VectorXi& get_start_map() const { return _start_map; }
    const Eigen::VectorXi& get_event_map() const { return _event_map; }
    const Eigen::VectorXd& get_scaling() const { return _scaling; }
    const Eigen::VectorXi& get_status() const { return _status; }
    const Eigen::VectorXd& get_event() const { return _event; }
    const Eigen::VectorXd& get_start() const { return _start; }

private:
    // Preprocessing results
    Eigen::VectorXd _start, _event, _scaling;
    Eigen::VectorXd raw_start;    // start times in native order
    Eigen::VectorXd raw_event;    // event times in native order
    Eigen::VectorXi raw_status;   // status in native order
    Eigen::VectorXi _status, _first, _last, _start_map, _event_map;
    Eigen::VectorXi event_order, start_order;
    
    // Buffers
    Eigen::VectorXd T_1_term, T_2_term, forward_scratch_buffer, w_avg_buffer, exp_w_buffer;
    Eigen::VectorXd grad_buffer, diag_hessian_buffer, diag_part_buffer, hess_matvec_buffer;
    Eigen::VectorXd linear_predictor, sample_weight;
    
    std::vector<Eigen::VectorXd> event_reorder_buffers;
    std::vector<Eigen::VectorXd> forward_cumsum_buffers;
    std::vector<Eigen::VectorXd> reverse_cumsum_buffers;
    std::vector<Eigen::VectorXd> risk_sum_buffers;

    double loglik_sat;
    bool have_start_times;
    bool _efron;

    void setup_buffers(int n);

    CoxContext get_context() const {
        return {event_order, start_order, _status, _first, _last, _scaling, _event_map, _start_map, have_start_times, _efron};
    }
};

class StratifiedCoxDeviance {
public:
    StratifiedCoxDeviance(const Eigen::VectorXd& start,
                          const Eigen::VectorXd& event,
                          const Eigen::VectorXi& status,
                          const Eigen::VectorXi& strata,
                          const Eigen::VectorXd& weight,
                          bool efron);

    double compute_deviance(const Eigen::VectorXd& eta,
                            const Eigen::VectorXd& sample_weight);

    void compute_hessian_matvec(const Eigen::VectorXd& arg,
                                Eigen::VectorXd& out);

    // As for CoxDeviance, stratum by stratum (each with its own seed derived
    // from seed).
    Eigen::MatrixXd sample_score(const std::string& method,
                                   int n_draws,
                                   std::uint64_t seed);
    Eigen::VectorXd wild_score(const Eigen::VectorXd& multipliers);
    Eigen::VectorXd score_residuals(const Eigen::VectorXd& x) const;
    // within each stratum, labelled by its value of strata
    BaselineHazard baseline_hazard() const;
    // comparing only pairs within a stratum: the sums of the strata's
    // numerators and denominators
    Concordance concordance(const Eigen::VectorXd& pred,
                            const Eigen::VectorXd& weight) const;
    Eigen::VectorXd schoenfeld_residuals(const Eigen::VectorXd& x) const;
    Eigen::VectorXd risk_set_variance(const Eigen::VectorXd& z) const;
    Eigen::VectorXd expected_events() const;

    const Eigen::VectorXd& get_gradient() const { return grad_buffer; }
    const Eigen::VectorXd& get_diag_hessian() const { return diag_hessian_buffer; }
    const Eigen::VectorXd& get_linear_predictor() const { return linear_predictor; }
    const Eigen::VectorXd& get_sample_weight() const { return sample_weight; }
    double get_loglik_sat() const { return loglik_sat; }

    const std::vector<int>& get_unique_strata() const { return unique_strata; }
    const std::vector<std::vector<int>>& get_stratum_indices() const { return stratum_indices; }
    const std::vector<std::shared_ptr<CoxDeviance>>& get_cox_devs() const { return cox_devs; }

private:
    std::vector<int> unique_strata;
    std::vector<std::vector<int>> stratum_indices;
    std::vector<std::shared_ptr<CoxDeviance>> cox_devs;

    Eigen::VectorXd grad_buffer;
    Eigen::VectorXd diag_hessian_buffer;
    Eigen::VectorXd linear_predictor;
    Eigen::VectorXd sample_weight;
    double loglik_sat;
};

// Harrell's concordance index (Harrell et al., 1982, JAMA 247, 2543-2546;
// Harrell, Lee and Mark, 1996, Statistics in Medicine 15, 361-387) of a
// predictor pred for (start, event] data, all n-vectors in the same order.
// For each failure i at time t, every j at risk at t (start_j < t <= event_j)
// is compared with i, except the other failures at t: the pair counts with
// weight w_i w_j, as concordant (1) if pred_i > pred_j and half if they are
// tied. This is R's survival::concordance(y ~ -pred, weights = w) and
// glmnet::Cindex. O(n log n), by the algorithm of the survival package's
// concordance vignette (Therneau and Atkinson): subjects are taken in
// decreasing time, with the weights of the risk set kept in a Fenwick
// (binary indexed) tree over the ranks of pred (Fenwick, 1994, Software:
// Practice and Experience 24, 327-336), counting pairs in order as Knight
// (1966, JASA 61, 436-439) does for Kendall's tau.
Concordance concordance(const Eigen::VectorXd& pred,
                        const Eigen::VectorXd& start,
                        const Eigen::VectorXd& event,
                        const Eigen::VectorXi& status,
                        const Eigen::VectorXd& weight);

// Standalone function declarations
void forward_cumsum(const Eigen::Ref<const Eigen::VectorXd> sequence,
		    Eigen::Ref<Eigen::VectorXd> output);

void reverse_cumsums(const Eigen::Ref<const Eigen::VectorXd> sequence,
                     Eigen::Ref<Eigen::VectorXd> event_buffer,
                     Eigen::Ref<Eigen::VectorXd> start_buffer,
                     const Eigen::Ref<const Eigen::VectorXi> event_order,
                     const Eigen::Ref<const Eigen::VectorXi> start_order,
		     bool do_event = false,
		     bool do_start = false);

void to_native_from_event(Eigen::Ref<Eigen::VectorXd> arg,
			  const Eigen::Ref<const Eigen::VectorXi> event_order,
			  Eigen::Ref<Eigen::VectorXd> reorder_buffer);

void to_event_from_native(const Eigen::Ref<const Eigen::VectorXd> arg,
                          const Eigen::Ref<const Eigen::VectorXi> event_order,
                          Eigen::Ref<Eigen::VectorXd> reorder_buffer);

void forward_prework(const Eigen::Ref<const Eigen::VectorXi> status,
                     const Eigen::Ref<const Eigen::VectorXd> w_avg,
                     const Eigen::Ref<const Eigen::VectorXd> scaling,
                     const Eigen::Ref<const Eigen::VectorXd> risk_sums,
                     int i,
                     int j,
                     Eigen::Ref<Eigen::VectorXd> moment_buffer,
		     const Eigen::Ref<const Eigen::VectorXd> arg,		     
                     bool use_w_avg = true);

double compute_sat_loglik(const Eigen::Ref<const Eigen::VectorXi> first,
			  const Eigen::Ref<const Eigen::VectorXi> last,
			  const Eigen::Ref<const Eigen::VectorXd> weight,
			  const Eigen::Ref<const Eigen::VectorXi> event_order,
			  const Eigen::Ref<const Eigen::VectorXi> status,
			  const Eigen::Ref<const Eigen::VectorXd> scaling,
			  Eigen::Ref<Eigen::VectorXd> W_status,
			  bool efron);

void sum_over_events(const CoxContext& ctx,
                     Eigen::Ref<Eigen::VectorXd> C_arg,
                     Eigen::Ref<Eigen::VectorXd> C_arg_scale,
                     Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
                     Eigen::Ref<Eigen::VectorXd> value_buffer);

void sum_over_risk_set(const Eigen::Ref<const Eigen::VectorXd> arg,
                       const CoxContext& ctx,
                       Eigen::Ref<Eigen::VectorXd> risk_sum_buffer,
                       Eigen::Ref<Eigen::VectorXd> event_cumsum,
                       Eigen::Ref<Eigen::VectorXd> start_cumsum);

double cox_dev(const Eigen::Ref<const Eigen::VectorXd> eta,
	       const Eigen::Ref<const Eigen::VectorXd> sample_weight,
	       const Eigen::Ref<const Eigen::VectorXd> exp_w,
               const CoxContext& ctx,
	       double loglik_sat,
	       Eigen::Ref<Eigen::VectorXd> T_1_term,
	       Eigen::Ref<Eigen::VectorXd> T_2_term,
	       Eigen::Ref<Eigen::VectorXd> grad_buffer,
	       Eigen::Ref<Eigen::VectorXd> diag_hessian_buffer,
	       Eigen::Ref<Eigen::VectorXd> diag_part_buffer,
	       Eigen::Ref<Eigen::VectorXd> w_avg_buffer,
               Eigen::Ref<Eigen::VectorXd> eta_event,
               Eigen::Ref<Eigen::VectorXd> w_event,
               Eigen::Ref<Eigen::VectorXd> exp_eta_w_event,
               Eigen::Ref<Eigen::VectorXd> risk_sums,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum0,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum1,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum2,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum3,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum4,
	       Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
               Eigen::Ref<Eigen::VectorXd> event_cumsum,
               Eigen::Ref<Eigen::VectorXd> start_cumsum);

void hessian_matvec(const Eigen::Ref<const Eigen::VectorXd> arg,
                    const CoxContext& ctx,
                    Eigen::Ref<Eigen::VectorXd> risk_sums,
                    Eigen::Ref<const Eigen::VectorXd> diag_part,
                    Eigen::Ref<const Eigen::VectorXd> w_avg,
                    Eigen::Ref<const Eigen::VectorXd> exp_w,
                    Eigen::Ref<Eigen::VectorXd> risk_sum_arg,
                    Eigen::Ref<Eigen::VectorXd> forward_cumsum0,
                    Eigen::Ref<Eigen::VectorXd> forward_cumsum1,
                    Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
                    Eigen::Ref<Eigen::VectorXd> event_cumsum,
                    Eigen::Ref<Eigen::VectorXd> start_cumsum,
                    Eigen::Ref<Eigen::VectorXd> hess_matvec_buffer);

void setup_preprocess(const Eigen::Ref<const Eigen::VectorXd> start,
                      const Eigen::Ref<const Eigen::VectorXd> event,
                      const Eigen::Ref<const Eigen::VectorXi> status,
                      const Eigen::Ref<const Eigen::VectorXd> weight,
                      Eigen::VectorXd& _start,
                      Eigen::VectorXd& _event,
                      Eigen::VectorXi& _status,
                      Eigen::VectorXi& _first,
                      Eigen::VectorXi& _last,
                      Eigen::VectorXd& _scaling,
                      Eigen::VectorXi& _start_map,
                      Eigen::VectorXi& _event_map,
                      Eigen::VectorXi& event_order,
                      Eigen::VectorXi& start_order);
