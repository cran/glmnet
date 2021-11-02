#pragma once
#include <functional>
#include <type_traits>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/util/iterator/counting_iterator.hpp>
#include <glmnetpp_bits/util/iterator/one_to_zero_iterator.hpp>
#include <Eigen/Core>

namespace glmnetpp {

/*
 * Base class for internal implementation of Gaussian elastic-net point solver.
 * This contains all the common interface and members across all versions of gaussian:
 *      - dense gaussian cov
 *      - dense gaussian naive
 *      - sparse gaussian cov
 *      - sparse gaussian naive
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalGaussianBase
{
    using value_t = ValueType;
    using index_t = IndexType;
    using bool_t = BoolType;

    template <class IAType
            , class XVType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalGaussianBase(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju)
        : thr_(thr)
        , maxit_(maxit)
        , nx_(nx)
        , a_(xv.size())
        , mm_(xv.size())
        , nlp_(nlp)
        , ia_(ia.data(), ia.size())
        , xv_(xv.data(), xv.size())
        , vp_(vp.data(), vp.size())
        , cl_(cl.data(), cl.rows(), cl.cols())
        , ju_(util::init_bvec<bool_t>::eval(ju))
    {
        a_.setZero(); 
        mm_.setZero();
    }

    bool is_warm() const { return !jz_; }
    bool is_warm_ever() const { return iz_; }
    void set_warm() { jz_ = false; }
    void reset_warm() { jz_ = true; }
    void set_warm_ever() { iz_ = true; }
    index_t& passes() { return nlp_; }
    void coord_desc_reset() { dlx_ = 0.0; }
    auto active_begin() const { return util::one_to_zero_iterator<index_t>(ia_.data()); }
    auto active_end() const { return util::one_to_zero_iterator<index_t>(ia_.data() + nin_); }
    constexpr auto all_begin() const { return util::counting_iterator<index_t>(0); }
    auto all_end() const { return util::counting_iterator<index_t>(a_.size()); }
    bool has_converged() const { return dlx_ < thr_; }
    bool has_reached_max_passes() const { return nlp_ > maxit_; }
    bool has_reached_max_active() const { return nin_ > nx_; }
    bool is_active(index_t j) const { return mm_(j) != 0; }

    template <class Iter
            , class UpdatePolicy
            , class SkipPolicy>
    constexpr void coord_desc(
            Iter begin,
            Iter end,
            UpdatePolicy update_pol,
            SkipPolicy skip_pol) const
    {
        for (; begin != end; ++begin) {
            auto k = *begin; 
            if (skip_pol(k)) continue;
            update_pol(k);
        } 
    }

    auto beta(index_t k) const { return a_(k); }

    void update_beta(index_t k, value_t ab, value_t dem, value_t gk) {
        auto ak = a_(k); 
        auto u = gk + ak * xv_(k); 
        auto v = std::abs(u) - vp_(k) * ab; 
        a_(k) = 0.0;
        if (v > 0.0) {
            a_(k) = std::max(cl_(0,k),
                    std::min(cl_(1,k),
                        std::copysign(v,u) / (xv_(k) + vp_(k) * dem)) );
        }
    }

    void update_active(index_t k)
    {
        ++nin_;
        if (has_reached_max_active()) {
            throw util::max_active_reached_error();
        }
        mm_(k) = nin_; 
        ia_(nin_-1) = k+1;
    }

    void update_rsq(index_t k, value_t beta_diff, value_t gk) { 
        rsq_ += beta_diff * (2.0 * gk - beta_diff * xv_(k));
    }

    void update_dlx(index_t k, value_t beta_diff) {
        dlx_ = std::max(xv_(k) * beta_diff * beta_diff, dlx_);
    }

    auto rsq() const { return rsq_; }
    auto n_active() const { return nin_; }

protected:
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;
    using ivec_t = Eigen::Matrix<index_t, Eigen::Dynamic, 1>;
    using bvec_t = util::bvec_t<bool_t>;
    using ju_t = typename std::conditional<
                    std::is_same<bool_t, bool>::value,
                    const bvec_t&,
                    Eigen::Map<const bvec_t> >::type;

    // internal non-captures
    bool iz_ = false;           // true if a partial fit was done with a previous lambda (warm ever)
    bool jz_ = true;            // true if requires a partial fit for current lambda (not warm)
    value_t dlx_ = 0.0;         // measures convergence of each coord-desc
    value_t rsq_ = 0.0;         // R^2
    const value_t thr_;         // threshold for convergence
    const index_t maxit_;       // max number of passes
    index_t nin_ = 0;           // number of active variables
    const index_t nx_;          // max number of active variables
    vec_t a_;                   // uncompressed beta
    ivec_t mm_;                 // index k is j if feature k is the jth feature active

    // captures
    index_t& nlp_;                      // number of passes
    Eigen::Map<ivec_t> ia_;             // active set (important that it's 1-indexed!)
    Eigen::Map<const vec_t> xv_;        // variance of columns of x
    Eigen::Map<const vec_t> vp_;        // penalties on features
    Eigen::Map<const mat_t> cl_;        // limits on each feature (2 x nvars)
    ju_t ju_;                           // exclusion type
};


/*
 * Base class for internal implementation of Gaussian covariance method.
 * This contains all the common interface and members for gaussian covariance methods:
 *      - dense gaussian cov
 *      - sparse gaussian cov
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalGaussianCovBase
        : ElnetPointInternalGaussianBase<
            ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianBase<
            ValueType, IndexType, BoolType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;

    template <class IAType
            , class GType
            , class XVType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalGaussianCovBase(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            GType& g,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju)
        : base_t(thr, maxit, nx, nlp, ia, xv, vp, cl, ju)
        , da_(g.size())
        , g_(g.data(), g.size())
        , c_(g.size(), nx)
    {}

    bool is_excluded(index_t j) const { return ju_[j] == 0; }

    constexpr void initialize(value_t, value_t, value_t) const { return; }

    template <class InitialFitIntType>
    constexpr bool initial_fit(InitialFitIntType f) const 
    { 
        bool converged = false, kkt_passed = false;
        std::tie(converged, kkt_passed) = f(); 
        return kkt_passed;
    }

    void update_beta(index_t k, value_t ab, value_t dem) {
        base_t::update_beta(k, ab, dem, g_(k));
    }

    void compress_active() {
        for (index_t j = 0; j < nin_; ++j) {
            da_(j) = a_(ia_(j)-1);
        }
    }

    void update_compressed_active() { 
        for (index_t j = 0; j < nin_; ++j) {
            da_(j) -= a_(ia_(j)-1);
        }
    }

    void update_grad_compressed_active() {
        for (index_t j = 0; j < a_.size(); ++j) {
            if (this->is_active(j)) continue;
            if (!is_excluded(j)) {
                g_(j) += da_.head(nin_).dot(c_.row(j).head(nin_));
            }
        }
    }

    void update_rsq(index_t k, value_t beta_diff) {
        base_t::update_rsq(k, beta_diff, g_(k));
    }

    void update_grad(index_t j, index_t k, value_t beta_diff) {
        g_(j) -= c_(j, mm_(k)-1) * beta_diff;
    }

    constexpr bool check_kkt(value_t) const { return true; }

protected:
    using typename base_t::vec_t;
    using typename base_t::mat_t;
    using base_t::ju_;
    using base_t::nin_;
    using base_t::a_;
    using base_t::ia_;
    using base_t::mm_;
    using base_t::xv_;

    vec_t da_;                  // compressed beta
    Eigen::Map<vec_t> g_;       // gradient
    mat_t c_;                   // X^TX cache but of dimension (nvars, max_nvars)
};

/*
 * Base class for internal implementation of Gaussian naive method.
 * This contains all the common interface and members for gaussian naive methods:
 *      - dense gaussian naive
 *      - sparse gaussian naive
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalGaussianNaiveBase
        : ElnetPointInternalGaussianBase<
            ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianBase<
            ValueType, IndexType, BoolType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;

    template <class IAType
            , class XVType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalGaussianNaiveBase(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju)
        : base_t(thr, maxit, nx, nlp, ia, xv, vp, cl, ju)
        , g_(ju.size())
        , ix_(ju.size(), false)
    {
        g_.setZero();
    }

    bool is_excluded(index_t j) const { return !ix_[j]; }

    void initialize(value_t beta, value_t alm, value_t alm0) { 
        auto tlam = beta * (2.0 * alm - alm0);
        for (index_t k = 0; k < g_.size(); ++k) {
            if (!is_excluded(k) || ju_[k] == 0) continue;
            if (g_(k) > tlam * vp_(k)) ix_[k] = true;
        }
    }

    template <class InitialFitIntType>
    constexpr bool initial_fit(InitialFitIntType f) const 
    { 
        // Keep doing initial fit until either doesn't converge or 
        // converged and kkt passed.
        while (1) {
            if (this->has_reached_max_passes()) { 
                throw util::maxit_reached_error();
            }
            bool converged = false, kkt_passed = false;
            std::tie(converged, kkt_passed) = f();
            if (!converged) break;
            if (kkt_passed) return true;
        }
        return false;
    }

    void update_beta(index_t k, value_t ab, value_t dem, value_t gk) {
        gk_cache_ = gk; 
        base_t::update_beta(k, ab, dem, gk_cache_);
    }

    void update_rsq(index_t k, value_t beta_diff) { 
        base_t::update_rsq(k, beta_diff, gk_cache_);
    }

    template <class GradFType>
    bool check_kkt(value_t ab, GradFType grad_f) {
        bool kkt_passed = true;
        for (index_t k = 0; k < g_.size(); ++k) {
            if (!is_excluded(k) || ju_[k] == 0) continue;
            g_(k) = std::abs(grad_f(k));
            if (g_(k) > ab * vp_(k)) { 
                ix_[k] = true; 
                kkt_passed = false;
            }
        }
        return kkt_passed;
    }

    const auto& abs_grad() const { return g_; }

protected:

    template <class GradFType>
    void compute_abs_grad(GradFType grad_f) {
        for (index_t j = 0; j < g_.size(); ++j) {
            if (ju_[j] == 0) continue;
            g_(j) = std::abs(grad_f(j));
        }
    }

    using typename base_t::vec_t;
    using typename base_t::mat_t;
    using base_t::ju_;
    using base_t::vp_;

    value_t gk_cache_ = 0.0;    // caches gradient at k when updating beta
    vec_t g_;                   // cwise-absolute gradient
    std::vector<bool> ix_;      // strong set indicators
};

} // namespace glmnetpp
