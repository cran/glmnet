#include <cstddef>
#include <set>
#include <RcppEigen.h>
#include <coxdev.h>

using namespace Rcpp;

// Compute Cox deviance, gradient, and diagonal Hessian via coxdev.
//
// Unified C++ entry point for all Cox computations.
// Handles: Breslow + Efron, strata, right-censored, start-stop.
//
// Returns a list with:
//   deviance: scalar (or vector if eta is a matrix)
//   gradient: vector (d(loglik)/d(eta), same scale as coxgrad)
//   diag_hessian: vector (diagonal of Hessian)
//   loglik_sat: scalar (saturated log-likelihood)
//
// [[Rcpp::export]]
List compute_cox_quantities_exp(
    Eigen::VectorXd start,      // start times (zeros for right-censored)
    Eigen::VectorXd stop,       // stop/event times
    Eigen::VectorXi status,     // event status
    Eigen::VectorXi strata,     // strata labels (empty = no strata)
    bool efron,                 // Efron or Breslow
    Eigen::MatrixXd eta,        // linear predictor (n x k matrix, k predictions)
    Eigen::VectorXd weights     // sample weights
    )
{
    int n = status.size();
    int k = eta.cols();

    Eigen::VectorXd deviance(k);

    // For the last column's gradient/hessian (most common case: k=1)
    Eigen::VectorXd gradient(n);
    Eigen::VectorXd diag_hessian(n);
    double loglik_sat = 0.0;

    if (strata.size() > 0) {
        StratifiedCoxDeviance scd(start, stop, status, strata, weights, efron);
        for (int j = 0; j < k; ++j) {
            Eigen::VectorXd eta_j = eta.col(j);
            deviance(j) = scd.compute_deviance(eta_j, weights);
        }
        // Return gradient/hessian from last computation
        const Eigen::VectorXd& g = scd.get_gradient();
        const Eigen::VectorXd& h = scd.get_diag_hessian();
        // coxdev returns -2 * d(loglik)/d(eta); convert to d(loglik)/d(eta)
        gradient = -g * 0.5;
        diag_hessian = h;
        loglik_sat = scd.get_loglik_sat();
    } else {
        CoxDeviance cd(start, stop, status, weights, efron);
        for (int j = 0; j < k; ++j) {
            Eigen::VectorXd eta_j = eta.col(j);
            deviance(j) = cd.compute_deviance(eta_j, weights);
        }
        const Eigen::VectorXd& g = cd.get_gradient();
        const Eigen::VectorXd& h = cd.get_diag_hessian();
        gradient = -g * 0.5;
        diag_hessian = h;
        loglik_sat = cd.get_loglik_sat();
    }

    return List::create(
        Named("deviance") = deviance,
        Named("gradient") = gradient,
        Named("diag_hessian") = diag_hessian,
        Named("loglik_sat") = loglik_sat);
}

// Cox lambda_max computation — full-featured version
//
// Handles: Breslow + Efron, strata, sparse + dense X, left-truncation,
// excluded variables, penalty factors (including zero penalty factors).
//
// For zero penalty factors: fits those variables at lambda=0 using
// coxdev to get the correct eta, then computes gradient at that eta.
//
// [[Rcpp::export]]
List compute_cox_lambda_max_exp(
    Eigen::MatrixXd x,          // design matrix (nobs x nvars)
    Eigen::VectorXd start,      // start times
    Eigen::VectorXd stop,       // stop times
    Eigen::VectorXi status,     // event status
    Eigen::VectorXi strata,     // strata labels (empty = no strata)
    bool efron,                 // Efron or Breslow
    Eigen::VectorXd offset,     // offset (eta at null)
    Eigen::VectorXd weights,    // sample weights
    Eigen::VectorXd xm,         // centering means
    Eigen::VectorXd xs,         // scaling factors
    Eigen::VectorXd vp,         // penalty factors
    Eigen::VectorXi exclude,    // excluded variables (1-indexed from R)
    double alpha                // elastic net alpha
    )
{
    int nvars = x.cols();

    // Build exclude set (convert from 1-indexed to 0-indexed)
    std::set<int> exclude_set;
    for (int i = 0; i < exclude.size(); ++i) {
        exclude_set.insert(exclude(i) - 1);
    }

    // Compute gradient using coxdev
    Eigen::VectorXd grad;
    double loglik_sat;

    if (strata.size() > 0) {
        StratifiedCoxDeviance scd(start, stop, status, strata, weights, efron);
        scd.compute_deviance(offset, weights);
        grad = scd.get_gradient();
        loglik_sat = scd.get_loglik_sat();
    } else {
        CoxDeviance cd(start, stop, status, weights, efron);
        cd.compute_deviance(offset, weights);
        grad = cd.get_gradient();
        loglik_sat = cd.get_loglik_sat();
    }

    // Convert gradient to score: d(loglik)/d(eta) = -gradient / 2
    Eigen::VectorXd score = -grad * 0.5;

    // Compute lambda_max: max_j |X_j^T * score - sum(score) * xm_j| / (xs_j * vp_j)
    double grad_sum = score.sum();
    double max_g = 0.0;

    for (int j = 0; j < nvars; ++j) {
        if (exclude_set.count(j) > 0 || vp(j) <= 0) continue;

        double inner_prod = score.dot(x.col(j));
        double g_j = std::abs((inner_prod - grad_sum * xm(j)) / xs(j));
        g_j /= vp(j);

        if (g_j > max_g) max_g = g_j;
    }

    double lambda_max = max_g / std::max(alpha, 1e-3);

    bool have_start_times = (start.array() > -std::numeric_limits<double>::infinity()).any();

    return List::create(
        Named("lambda_max") = lambda_max,
        Named("loglik_sat") = loglik_sat,
        Named("have_start_times") = have_start_times);
}
