// Tests of concordance (and CoxDeviance::concordance, with strata) against
// the O(n^2) count from the definition of Harrell's C-index, on data with
// tied times, tied predictions, start times and zero weights.
//
// Build and run from the top of the repository:
//   c++ -std=c++17 -O2 -Ieigen -Icpp/include cpp/src/coxdev.cpp \
//       cpp/tests/test_concordance.cpp -o test_concordance
//   ./test_concordance

#include "coxdev.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <string>

namespace {

int failures = 0;

void check(double error, double tol, const std::string& what) {
  if (!(error <= tol)) {
    std::printf("FAIL %s: error %.3g > %.3g\n", what.c_str(), error, tol);
    ++failures;
  }
}

// for each failure i at t and each j at risk at t (start_j < t <= event_j)
// that is not a failure at t, in the same stratum: w_i w_j counted, and
// w_i w_j (1 + sign(pred_i - pred_j)) / 2 concordant
Concordance brute_force(const Eigen::VectorXd& pred, const Eigen::VectorXd& start,
                        const Eigen::VectorXd& event, const Eigen::VectorXi& status,
                        const Eigen::VectorXd& weight, const Eigen::VectorXi& strata) {
  Concordance c{0.0, 0.0};
  int n = pred.size();
  for (int i = 0; i < n; ++i) {
    if (!status(i)) continue;
    double t = event(i);
    for (int j = 0; j < n; ++j) {
      if (strata(j) != strata(i) || !(start(j) < t && t <= event(j))) continue;
      if (status(j) && event(j) == t) continue;
      double w = weight(i) * weight(j);
      double sign = (pred(i) > pred(j)) - (pred(i) < pred(j));
      c.numerator += w * (1 + sign) / 2;
      c.denominator += w;
    }
  }
  return c;
}

}  // namespace

int main() {
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::normal_distribution<double> norm(0.0, 1.0);
  int cases = 0;
  for (int rep = 0; rep < 40; ++rep) {
    int n = 2 + rep * 3;
    Eigen::VectorXd pred(n), start(n), event(n), weight(n);
    Eigen::VectorXi status(n), strata(n), no_strata = Eigen::VectorXi::Zero(n);
    for (int i = 0; i < n; ++i) {
      event(i) = 1 + std::floor(unif(rng) * 9);              // ties in time
      status(i) = unif(rng) < 0.6;
      start(i) = std::min(std::floor(event(i) * unif(rng)), event(i) - 1);
      weight(i) = unif(rng) < 0.1 ? 0.0 : 2 * unif(rng);      // some zero
      pred(i) = std::round(norm(rng) * 10) / 10;              // ties in pred
      strata(i) = i % 3;
    }
    Eigen::VectorXd no_start = Eigen::VectorXd::Constant(n, -std::numeric_limits<double>::infinity());
    for (int with_start = 0; with_start < 2; ++with_start) {
      const Eigen::VectorXd& s = with_start ? start : no_start;
      std::string name = "rep " + std::to_string(rep) + (with_start ? " start" : "");

      Concordance c = concordance(pred, s, event, status, weight);
      Concordance b = brute_force(pred, s, event, status, weight, no_strata);
      check(std::abs(c.numerator - b.numerator), 1e-10 * (1 + b.numerator), name + ": numerator");
      check(std::abs(c.denominator - b.denominator), 1e-10 * (1 + b.denominator), name + ": denominator");

      CoxDeviance cox(s, event, status, weight, true);
      Concordance cc = cox.concordance(pred, weight);
      check(std::abs(cc.numerator - b.numerator) + std::abs(cc.denominator - b.denominator),
            1e-10 * (1 + b.denominator), name + ": CoxDeviance");

      StratifiedCoxDeviance scox(s, event, status, strata, weight, true);
      Concordance cs = scox.concordance(pred, weight);
      Concordance bs = brute_force(pred, s, event, status, weight, strata);
      check(std::abs(cs.numerator - bs.numerator) + std::abs(cs.denominator - bs.denominator),
            1e-10 * (1 + bs.denominator), name + ": strata");

      // -pred reverses the concordance of the untied pairs: C(pred) + C(-pred) = 1
      Concordance neg = concordance(-pred, s, event, status, weight);
      if (c.denominator > 0) {
        check(std::abs(c.value() + neg.value() - 1), 1e-12, name + ": C(pred) + C(-pred)");
      }
      ++cases;
    }
  }
  // a constant predictor ties every pair: 1/2
  Eigen::VectorXd e(4), w = Eigen::VectorXd::Ones(4), p = Eigen::VectorXd::Zero(4);
  Eigen::VectorXi d(4);
  e << 1, 2, 3, 4;
  d << 1, 1, 0, 1;
  Eigen::VectorXd s = Eigen::VectorXd::Constant(4, -std::numeric_limits<double>::infinity());
  check(std::abs(concordance(p, s, e, d, w).value() - 0.5), 1e-15, "constant predictor");
  std::printf("%d cases, %d failures\n", cases, failures);
  return failures == 0 ? 0 : 1;
}
