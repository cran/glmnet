#pragma once
#include <functional>
#include <type_traits>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_point/internal/base.hpp>

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
    : ElnetPointInternalBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::vec_t;
    using typename base_t::ivec_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

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
            const JUType& ju,
            value_t rsq = 0.0)
        : base_t(thr, maxit, nx, nlp, ia, vp, cl, ju)
        , rsq_(rsq)
        , xv_(xv.data(), xv.size())
    {}

    GLMNETPP_STRONG_INLINE bool is_warm_ever() const { return iz_; }
    GLMNETPP_STRONG_INLINE void set_warm_ever() { iz_ = true; }

    GLMNETPP_STRONG_INLINE
    void update_dlx(index_t k, value_t beta_diff) {
        base_t::update_dlx(beta_diff, xv_(k));
    }

    GLMNETPP_STRONG_INLINE constexpr void update_intercept() const {}

    GLMNETPP_STRONG_INLINE auto rsq() const { return rsq_; }
    GLMNETPP_STRONG_INLINE auto rsq_prev() const { return rsq_prev_; }

    /* Static interface */

    GLMNETPP_STRONG_INLINE
    static void
    update_rsq(value_t& rsq, value_t beta_diff, value_t gk, value_t x_var) { 
        rsq += beta_diff * (2.0 * gk - beta_diff * x_var);
    }

protected:
    using base_t::update_dlx;
    using base_t::update_intercept;

    GLMNETPP_STRONG_INLINE auto& rsq() { return rsq_; }
    GLMNETPP_STRONG_INLINE void initialize() { rsq_prev_ = rsq_; }

    GLMNETPP_STRONG_INLINE
    void update_rsq(index_t k, value_t beta_diff, value_t gk) { 
        update_rsq(rsq_, beta_diff, gk, xv_(k));
    }

    GLMNETPP_STRONG_INLINE
    auto x_var(index_t i) const { return xv_[i]; }

private:
    // internal non-captures
    bool iz_ = false;           // true if a partial fit was done with a previous lambda (warm ever)
    value_t rsq_;               // R^2
    value_t rsq_prev_ = 0.0;    // previous R^2

    // captures
    Eigen::Map<const vec_t> xv_;        // variance of columns of x
};

/*
 * Base class for internal implementation of Gaussian univariate-response methods.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalGaussianUniBase
    : ElnetPointInternalGaussianBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::vec_t;
    using typename base_t::value_t;
    using typename base_t::index_t;

    template <class IAType
            , class XVType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalGaussianUniBase(
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
        , a_(xv.size())
    {
        a_.setZero(); 
    }

    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t ab, value_t dem, value_t gk) {
        const auto& cl = this->endpts();
        base_t::update_beta(
                a_(k), gk, this->x_var(k), this->penalty()(k),
                cl(0,k), cl(1,k), ab, dem);
    }

public:
    GLMNETPP_STRONG_INLINE auto beta(index_t k) const { return a_(k); }

private:
    vec_t a_;                   // uncompressed beta
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
        : ElnetPointInternalGaussianUniBase<
            ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianUniBase<
            ValueType, IndexType, BoolType>;
protected:
    using typename base_t::vec_t;
    using typename base_t::mat_t;

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

    GLMNETPP_STRONG_INLINE bool is_excluded(index_t j) const { return !this->exclusion()[j]; }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PointPackType&) { base_t::initialize(); }

    template <class InitialFitIntType>
    GLMNETPP_STRONG_INLINE
    constexpr bool initial_fit(InitialFitIntType f) const 
    { 
        bool converged = false, kkt_passed = false;
        std::tie(converged, kkt_passed) = f(); 
        return kkt_passed;
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        base_t::update_beta(k, pack.ab, pack.dem, g_(k));
    }

    GLMNETPP_STRONG_INLINE
    void compress_active() {
        index_t j = 0;
        std::for_each(this->active_begin(), this->active_end(),
                [&](auto k) { da_(j++) = this->beta(k); });
    }

    GLMNETPP_STRONG_INLINE
    void update_compressed_active() { 
        index_t j = 0;
        std::for_each(this->active_begin(), this->active_end(),
                [&](auto k) { da_(j++) -= this->beta(k); });
    }

    GLMNETPP_STRONG_INLINE
    void update_grad_compressed_active() {
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                [&](auto j) {
                    auto n_act = this->n_active();
                    g_(j) += da_.head(n_act).dot(c_.row(j).head(n_act));
                },
                [&](auto j) { return this->is_active(j) || this->is_excluded(j); });
    }

    GLMNETPP_STRONG_INLINE
    void update_rsq(index_t k, value_t beta_diff) {
        base_t::update_rsq(k, beta_diff, g_(k));
    }

    GLMNETPP_STRONG_INLINE
    void update_grad(index_t j, index_t k, value_t beta_diff) {
        g_(j) -= c_(j, this->active_idx(k)) * beta_diff;
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    constexpr bool check_kkt(const PointPackType&) const { return true; }

protected:
    template <class XTXFType>
    GLMNETPP_STRONG_INLINE
    void update_active(index_t k, XTXFType xtx_f) {
        base_t::update_active(k);
        for (auto it = this->all_begin(); it != this->all_end(); ++it) {
            auto j = *it;
            if (this->is_excluded(j)) continue;

            // Note: j == k case is after the is_active(j) check in legacy code.
            // Since base_t::update_active adds mm_(k),
            // it will now be considered active,
            // so we have to first check this case.
            if (j == k) { 
                c_(j, this->n_active()-1) = this->x_var(j); 
                continue; 
            }
            if (this->is_active(j)) {
                c_(j, this->n_active()-1) = c_(k, this->active_idx(j)); 
                continue;
            }
            c_(j, this->n_active()-1) = xtx_f(j, k);
        }
    }

private:
    vec_t da_;                  // compressed beta
    Eigen::Map<vec_t> g_;       // gradient (not absolute gradient)
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
        : ElnetPointInternalGaussianUniBase<
            ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianUniBase<
            ValueType, IndexType, BoolType>;

protected:
    using typename base_t::vec_t;
    using typename base_t::mat_t;

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
    
    using base_t::update_intercept;

    GLMNETPP_STRONG_INLINE bool is_excluded(index_t j) const { return !ix_[j]; }

    template <class InitialFitIntType>
    GLMNETPP_STRONG_INLINE
    constexpr bool initial_fit(InitialFitIntType f) const {
        return initial_fit([&]() { return this->has_reached_max_passes(); }, f);
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PointPackType& pack) {
        base_t::initialize();
        initialize_strong_set(pack);
    }

    GLMNETPP_STRONG_INLINE
    void update_rsq(index_t k, value_t beta_diff) { 
        base_t::update_rsq(k, beta_diff, gk_cache_);
    }

    template <class AbsGradFType>
    GLMNETPP_STRONG_INLINE
    bool check_kkt(value_t ab, AbsGradFType abs_grad_f) {
        auto skip_f = [&](auto k) { return !is_excluded(k) || !this->exclusion()[k]; };
        update_abs_grad(g_, abs_grad_f, skip_f);
        return check_kkt(g_, this->penalty(), ix_, ab, skip_f);
    }

    const auto& abs_grad() const { return g_; }

    /* Static interface */

    template <class RType, class XType>
    GLMNETPP_STRONG_INLINE
    static void
    update_resid(
            RType&& r,
            value_t beta_diff,
            const XType& x) 
    {
        r -= beta_diff * x;
    }

    template <class RType, class VType>
    GLMNETPP_STRONG_INLINE
    static value_t 
    update_intercept(
            value_t& intercept,
            RType&& r,
            value_t& dlx,
            bool intr,
            value_t r_sum,
            value_t var,
            const VType& v)
    {
        auto d = base_t::update_intercept(intercept, dlx, intr, r_sum, var);
        if (d) update_resid(r, d, v);
        return d;
    }

    template <class HasReachedMaxPassesType, class InitialFitIntType>
    GLMNETPP_STRONG_INLINE
    constexpr static bool initial_fit(
            HasReachedMaxPassesType has_reached_max_passes, 
            InitialFitIntType f) 
    { 
        // Keep doing initial fit until either doesn't converge or 
        // converged and kkt passed.
        while (1) {
            if (has_reached_max_passes()) { 
                throw util::maxit_reached_error();
            }
            bool converged = false, kkt_passed = false;
            std::tie(converged, kkt_passed) = f();
            if (!converged) break;
            if (kkt_passed) return true;
        }
        return false;
    }

    /*
     * Updates absolute gradient abs_grad by iterating through each element
     * and assigning compute_grad_f(k). Iteration skips over k whenever skip_f(k) is true.
     */
    template <class AbsGradType, class ComputeAbsGradFType, class SkipFType>
    GLMNETPP_STRONG_INLINE
    static void update_abs_grad(
            AbsGradType&& abs_grad,
            ComputeAbsGradFType compute_abs_grad_f,
            SkipFType skip_f) 
    {
        base_t::for_each_with_skip(
                util::counting_iterator<index_t>(0), 
                util::counting_iterator<index_t>(abs_grad.size()),
                [&](index_t j) { abs_grad(j) = compute_abs_grad_f(j); },
                skip_f);
    }

    /*
     * Checks KKT condition and computes strong map. See base_t::compute_strong_map;
     * Returns true if no update occured (KKT all passed).
     */
    template <class AbsGradType
            , class PenaltyType
            , class StrongMapType
            , class SkipFType>
    GLMNETPP_STRONG_INLINE
    static bool check_kkt(
            AbsGradType&& abs_grad,
            const PenaltyType& penalty,
            StrongMapType&& strong_map,
            value_t l1_regul,
            SkipFType skip_f) {
        return !base_t::compute_strong_map(abs_grad, penalty, strong_map, l1_regul, skip_f);
    }

    template <class AbsGradType
            , class PenaltyType
            , class StrongMapType
            , class FType
            , class SkipFType>
    GLMNETPP_STRONG_INLINE
    static bool check_kkt(
            AbsGradType&& abs_grad,
            const PenaltyType& penalty,
            StrongMapType&& strong_map,
            value_t l1_regul,
            FType f,
            SkipFType skip_f) {
        return !base_t::compute_strong_map(abs_grad, penalty, strong_map, l1_regul, f, skip_f);
    }

protected:
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t ab, value_t dem, value_t gk) {
        gk_cache_ = gk; 
        base_t::update_beta(k, ab, dem, gk_cache_);
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void initialize_strong_set(const PointPackType& pack) 
    {
        base_t::compute_strong_map(
                g_, this->penalty(), ix_, 
                pack.elastic_prop(), pack.lmda(), pack.prev_lmda(),
                [&](auto k) { return !is_excluded(k) || !this->exclusion()[k]; }
                );
    }

    template <class AbsGradFType>
    GLMNETPP_STRONG_INLINE
    void construct(AbsGradFType abs_grad_f) {
        update_abs_grad(g_, abs_grad_f,
                [&](auto j) { return !this->exclusion()[j]; });
    }

private:
    value_t gk_cache_ = 0.0;    // caches gradient at k when updating beta
    vec_t g_;                   // cwise-absolute gradient
    std::vector<bool> ix_;      // strong set indicators
};

/*
 * Base class for internal implementation of multi-response Gaussian elastic-net point solver.
 * This contains all the common interface and members across all versions of multi-response gaussian:
 *      - dense gaussian multi
 *      - sparse gaussian multi
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalGaussianMultiBase 
    : ElnetPointInternalGaussianBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;

    // hide base protected members
    using base_t::initialize;

protected:
    using typename base_t::vec_t;
    using typename base_t::mat_t;
    using typename base_t::ivec_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class XVType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternalGaussianMultiBase(
            value_t thr,
            index_t maxit,
            index_t nr,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            value_t ys0,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(thr * ys0 / nr, maxit, nx, nlp, ia, xv, vp, cl /* actually won't be used */, ju, ys0)
        , bnorm_thr_(int_param.bnorm_thr)
        , bnorm_mxit_(int_param.bnorm_mxit)
        , ys0_(ys0)
        , nr_(nr)
        , a_(nr, xv.size())
        , g_curr_(nr)
        , g_next_(nr)
        , b_diff_(nr)
        , g_(xv.size())
        , ix_(xv.size(), false)
        , isc_(nr)
        , cl_(cl.data(), cl.size())
    {
        a_.setZero(); 
    }

    GLMNETPP_STRONG_INLINE auto& beta_buffer() { return b_diff_; }
    GLMNETPP_STRONG_INLINE const auto& abs_grad() const { return g_; }
    GLMNETPP_STRONG_INLINE bool is_excluded(index_t k) const { return !ix_[k]; }
    GLMNETPP_STRONG_INLINE auto beta(index_t k) const { return a_.col(k); }
    GLMNETPP_STRONG_INLINE auto beta(index_t k, index_t l) const { return a_(k, l); }
    constexpr GLMNETPP_STRONG_INLINE void update_intercept() const {}

    template <class DiffType>
    GLMNETPP_STRONG_INLINE void update_dlx(index_t k, const DiffType& diff) {
        base_t::update_dlx(diff, this->x_var(k));
    }

    /*
     * Slightly different formula for updating residual.
     * Our rsq starts at ys0_ and decrements down.
     * Later the rsq that we eventually save will be 1-rsq_/ys0_.
     */
    template <class DiffType>
    void update_rsq(index_t k, const DiffType& diff) 
    {
        this->rsq() -= (diff.array() * (2.0 * g_curr_ - this->x_var(k) * diff).array()).sum();
    }

    template <class InitialFitIntType>
    GLMNETPP_STRONG_INLINE
    constexpr bool initial_fit(InitialFitIntType f) const { 
        return gaussian_naive_t::initial_fit([&]() { return this->has_reached_max_passes(); }, f);
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE 
    void initialize(const PointPackType& pack) {
        base_t::initialize();
        base_t::compute_strong_map(
                g_, this->penalty(), ix_, pack.elastic_prop(), pack.lmda(), pack.prev_lmda(),
                [&](index_t k) { return !is_excluded(k) || !this->exclusion()[k]; });
    }

    template <class GradFType>
    GLMNETPP_STRONG_INLINE
    void update_beta(
            index_t k, 
            value_t l1_regul, 
            value_t l2_regul, 
            GradFType grad_f)
    {
        Eigen::Map<const mat_t> cl_slice(
                cl_.data() + k * 2 * nr_, 2, nr_);
        update_beta(k, a_.col(k), this->x_var(k), this->penalty()(k),
                g_curr_, g_next_, l1_regul, l2_regul, bnorm_thr_, bnorm_mxit_,
                isc_, cl_slice, grad_f);
    }

    /* Static interface */

    template <class AType
            , class GCurrType
            , class GNextType
            , class ISCType
            , class CLType
            , class GradFType>
    GLMNETPP_STRONG_INLINE
    static void update_beta(
            index_t k, 
            AType&& ak,
            value_t xvk,
            value_t vpk,
            GCurrType&& g_curr,
            GNextType&& g_next,
            value_t l1_regul, 
            value_t l2_regul, 
            value_t bnorm_thr,
            index_t bnorm_mxit,
            ISCType&& isc,
            const CLType& cl,
            GradFType grad_f)
    {
        grad_f(k, g_curr); 
        g_next = g_curr + xvk * ak;
        auto gkn = g_next.norm();
        auto u = 1.0 - l1_regul * vpk / gkn;
        if (u <= 0.0) { ak.setZero(); }
        else {
            ak = (u/(xvk + l2_regul * vpk)) * g_next;
            chkbnds(g_next, gkn, xvk, cl,
                    l2_regul*vpk, l1_regul*vpk, ak, isc, bnorm_thr, bnorm_mxit);
        }
    }

protected:
    template <class AbsGradFType>
    GLMNETPP_STRONG_INLINE
    void construct(AbsGradFType abs_grad_f) {
        gaussian_naive_t::update_abs_grad(g_, 
                [&](index_t j) { return abs_grad_f(j, g_curr_); },
                [&](auto j) { return !this->exclusion()[j]; });
    }

    template <class AbsGradFType>
    GLMNETPP_STRONG_INLINE
    bool check_kkt(value_t ab, AbsGradFType abs_grad_f) {
        auto skip_f = [&](auto k) { return !is_excluded(k) || !this->exclusion()[k]; };
        gaussian_naive_t::update_abs_grad(g_, 
                [&](index_t j) { return abs_grad_f(j, g_curr_); }, skip_f);
        return gaussian_naive_t::check_kkt(g_, this->penalty(), ix_, ab, skip_f);
    }

private:
    /*
     * TODO: Document what this is doing.
     * This and the next routine are only needed for multi-response stuff.
     * Maybe move this in a multi-only static interface?
     */
    template <class GKType
            , class CLType
            , class AType
            , class ISCType>
    GLMNETPP_STRONG_INLINE
    static void chkbnds(
            const GKType& gk,
            value_t gkn,
            value_t xv,
            const CLType& cl,
            value_t al1,
            value_t al2,
            AType& a,
            ISCType&& isc,
            value_t thr,
            index_t mxit) 
    {
        auto al1p = 1.0 + al1/xv;
        auto al2p = al2/xv;
        isc.setZero();
        auto gsq = gkn*gkn;
        auto asq = a.squaredNorm();
        double usq = 0.0;
        double u = 0.0;
        int kn = -1;
        while (1) {
            double vmx = 0.0;
            for (int k = 0; k < a.size(); ++k) {
                auto v = std::max(a(k)-cl(1,k), cl(0,k)-a(k));
                if (v > vmx) { vmx = v; kn = k; }
            }
            if (vmx <= 0.0) break;
            if (isc(kn)) break;
            gsq -= gk(kn)*gk(kn);
            // numerical stability: take abs
            auto g = std::sqrt(std::abs(gsq))/xv;
            if (a(kn) < cl(0, kn)) u = cl(0, kn);
            if (a(kn) > cl(1, kn)) u = cl(1, kn);
            usq += u*u;
            double b = 0.0;
            if (usq == 0.0) b = std::max(0., (g-al2p)/al1p);
            else {
                // numerical stability: take abs
                auto b0 = std::sqrt(std::abs(asq - a(kn) * a(kn))); 
                b = bnorm(b0, al1p, al2p, g, usq, thr, mxit);
            }
            asq = usq + b*b;
            if (asq <= 0.0) { a.setZero(); break; }
            a(kn) = u;
            isc(kn) = 1;
            auto f = 1.0/(xv * (al1p+al2p/std::sqrt(asq)));
            for (int j = 0; j < a.size(); ++j) {
                if (!isc(j)) a(j) = f * gk(j);
            }
        }
    }

    /*
     * TODO: Document what this is doing.
     */
    GLMNETPP_STRONG_INLINE
    static value_t bnorm(
            value_t b0,
            value_t al1p,
            value_t al2p,
            value_t g,
            value_t usq,
            value_t thr,
            index_t mxit)
    {
        auto b = b0;
        auto zsq = b*b + usq;
        if (zsq <= 0.0) { return 0.0; }
        auto z = std::sqrt(zsq);
        auto f = b*(al1p+al2p/z)-g;
        int it = 0;
        for (; it < mxit; ++it) {
            b -= f/(al1p+al2p*usq/(z*zsq));
            zsq = b*b + usq;
            if (zsq <= 0.0) { return 0.0; }
            z = std::sqrt(zsq);
            f = b*(al1p+al2p/z)-g;
            if (std::abs(f) <= thr) break;
            if (b <= 0.0) { b = 0.0; break; }
        }
        if (it >= mxit) throw util::bnorm_maxit_reached_error();
        return b;
    }

    const value_t bnorm_thr_;
    const index_t bnorm_mxit_;
    const value_t ys0_; // frobenius norm of y_ after standardizing
    const index_t nr_;  // number of responses
    mat_t a_;       // matrix of coefficients
    vec_t g_curr_;  // g_curr_(k) = y^T x_col_k
    vec_t g_next_;  // g_next_(k) = y^T x_col_k + beta_diff_k * xv(k)
    vec_t b_diff_;  // temporarily stores beta difference during CD
    vec_t g_;       // absolute gradient
    std::vector<bool> ix_; // strong set
    ivec_t isc_;    // buffer for efficiency during chkbnds
    Eigen::Map<const vec_t> cl_;    // override base capture of cl
                                    // this class specifically expects cl to be 3-d Array, which we capture as a vector.
};

/*
 * Base class for internal implementation of Gaussian general WLS elastic-net point solver.
 * This contains all the common interface and members across all versions of the general WLS solver:
 *      - dense gaussian wls
 *      - sparse gaussian wls
 * This class is intended for the general programmable method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalGaussianWLSBase 
        : ElnetPointInternalBaseViewer<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBaseViewer<ValueType, IndexType, BoolType>;
    using gaussian_t = ElnetPointInternalGaussianBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::ivec_t;
    using typename base_t::vec_t;
    using typename base_t::mat_t;

public:
    template <class RType
            , class XVType
            , class VType
            , class JUType
            , class VPType
            , class CLType
            , class AType
            , class GType
            , class IAType
            , class IYType
            , class MMType>
    ElnetPointInternalGaussianWLSBase(
        value_t alm0,
        value_t almc,
        value_t alpha,
        RType& r,
        XVType& xv,
        const VType& v,
        bool intr,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        index_t nx,
        value_t thr,
        index_t maxit,
        AType& a,
        value_t& aint,
        GType& g,
        IAType& ia,
        IYType& iy,
        index_t& iz,
        MMType& mm,
        index_t& nino,
        value_t& rsqc,
        index_t& nlp)
        : base_t(thr, maxit, nino, nx, nlp, ia, vp, cl, ju)
        , lmda_(almc)
        , prev_lmda_(alm0)
        , alpha_(alpha)
        , l1_regul_(almc * alpha)
        , l2_regul_(almc * (1.0 - alpha))
        , xmz_(v.sum())
        , intr_(intr)
        , iz_(iz)
        , rsq_(rsqc)
        , r_(r.data(), r.size())
        , xv_(xv.data(), xv.size())
        , v_(v.data(), v.size())
        , a_(a.data(), a.size())
        , a0_(aint)
        , g_(g.data(), g.size())
        , ix_(iy.data(), iy.size())
    {
        base_t::construct(mm);
    }
    
    GLMNETPP_STRONG_INLINE bool is_warm_ever() const { return iz_; }
    GLMNETPP_STRONG_INLINE void set_warm_ever() const { iz_ = true; }
    GLMNETPP_STRONG_INLINE bool is_excluded(index_t j) const { return !strong_map()[j]; }
    GLMNETPP_STRONG_INLINE value_t beta(index_t k) const { return a_(k); }

    /*
     * Tries initial fit until either it did not converge or both converged and kkt conditions are all satisfied.
     * The only difference from naive/multi is that it does not do an initial check of whether max iter has reached.
     */
    template <class InitialFitIntType>
    GLMNETPP_STRONG_INLINE
    constexpr bool initial_fit(InitialFitIntType f) const { 
        return gaussian_naive_t::initial_fit([&]() { return this->has_reached_max_passes(); }, f);
    }

    GLMNETPP_STRONG_INLINE
    void update_dlx(index_t k, value_t diff) {
        base_t::update_dlx(diff, xv_(k));
    }

    GLMNETPP_STRONG_INLINE
    void update_rsq(index_t k, value_t diff) {
        gaussian_t::update_rsq(rsq_, diff, gk_, xv_(k));
    }

protected:
    GLMNETPP_STRONG_INLINE auto& resid() { return r_; }
    GLMNETPP_STRONG_INLINE const auto& resid() const { return r_; }
    GLMNETPP_STRONG_INLINE const auto& weight() const { return v_; }
    GLMNETPP_STRONG_INLINE auto new_weight_sum() const { return xmz_; }
    
    template <class ComputeXVFType>
    GLMNETPP_STRONG_INLINE
    void initialize(ComputeXVFType compute_xv_f) {
        base_t::compute_strong_map(
                abs_grad(), this->penalty(), strong_map(), alpha_, lmda_, prev_lmda_, 
                [&](index_t j) { xv_(j) = compute_xv_f(j); },
                [&](index_t j) { return !is_excluded(j) || !this->exclusion()[j]; });
    }

    template <class XVFType, class AbsGradFType>
    GLMNETPP_STRONG_INLINE
    bool check_kkt(XVFType xv_f, AbsGradFType abs_grad_f) {
        auto skip_f = [&](index_t k) { return !is_excluded(k) || !this->exclusion()[k]; };
        gaussian_naive_t::update_abs_grad(g_, abs_grad_f, skip_f);
        return gaussian_naive_t::check_kkt(g_, this->penalty(), ix_, l1_regul_, 
                [&](index_t j) { xv_(j) = xv_f(j); }, skip_f);
    }

    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t gk) {
        gk_ = gk;
        const auto& cl = this->endpts();
        base_t::update_beta(
                a_(k), gk_, xv_(k), this->penalty()(k), 
                cl(0,k), cl(1,k), l1_regul_, l2_regul_);
    }

    GLMNETPP_STRONG_INLINE
    auto update_intercept(value_t r_sum) {
        auto d = gaussian_naive_t::update_intercept(
                a0_, this->resid(), this->convg_measure(), intr_,
                r_sum, xmz_, v_);
        if (d) gaussian_t::update_rsq(rsq_, d, r_sum, xmz_);
        return d;
    }

    template <class ComputeXVFType, class AbsGradFType>
    GLMNETPP_STRONG_INLINE
    void construct(ComputeXVFType compute_xv_f,
                   AbsGradFType abs_grad_f) 
    {
        gaussian_naive_t::update_abs_grad(g_, abs_grad_f,
                [&](index_t j) { return !this->exclusion()[j]; });
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                [&](index_t j) { xv_(j) = compute_xv_f(j); },
                [&](index_t j) { return is_excluded(j); });
    }

private:
    GLMNETPP_STRONG_INLINE auto& abs_grad() const { return g_; }
    GLMNETPP_STRONG_INLINE auto& strong_map() { return ix_; }
    GLMNETPP_STRONG_INLINE const auto& strong_map() const { return ix_; }

    value_t gk_ = 0.0; // caches current gradient
    const value_t lmda_;
    const value_t prev_lmda_;
    const value_t alpha_;
    const value_t l1_regul_;
    const value_t l2_regul_;
    const value_t xmz_;
    const bool intr_;
    index_t& iz_;
    value_t& rsq_;
    Eigen::Map<vec_t> r_;       // scaled residual vector
    Eigen::Map<vec_t> xv_;      // weighted variance of x-cols but not computed yet.
                                // Compute on-the-fly as features enter strong set.
    Eigen::Map<const vec_t> v_; // square-root of the weights.
    Eigen::Map<vec_t> a_;       // coefficients
    value_t& a0_;               // intercept
    Eigen::Map<vec_t> g_;       // absolute gradient
    Eigen::Map<ivec_t> ix_;     // strong map
};

} // namespace glmnetpp
