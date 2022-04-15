#pragma once
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/util/iterator/counting_iterator.hpp>
#include <glmnetpp_bits/util/iterator/one_to_zero_iterator.hpp>
#include <glmnetpp_bits/util/type_traits.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace glmnetpp {

/*
 * Base class for all static interface
 */
template <class ValueType, class IndexType>
struct ElnetPointInternalStaticBase
{
protected:
    using value_t = ValueType;
    using index_t = IndexType;

    /* Helper routines to encapsulate the logic of some common routines in a generic way */

    template <class Iter
            , class UpdatePolicy
            , class SkipPolicy>
    GLMNETPP_STRONG_INLINE
    static constexpr void 
    for_each_with_skip(
            Iter begin,
            Iter end,
            UpdatePolicy update_pol,
            SkipPolicy skip_pol) 
    {
        for (; begin != end; ++begin) {
            auto k = *begin; 
            if (skip_pol(k)) continue;
            update_pol(k);
        } 
    }

    GLMNETPP_STRONG_INLINE
    static void 
    update_dlx(value_t& dlx, value_t beta_diff, value_t x_var) {
        dlx = std::max(x_var * beta_diff * beta_diff, dlx);
    }

    GLMNETPP_STRONG_INLINE
    static void 
    update_beta(
            value_t& a, 
            value_t gk, 
            value_t x_var,
            value_t penalty,
            value_t a_min,
            value_t a_max,
            value_t l1_regul, 
            value_t l2_regul) {
        auto a_copy = a; 
        auto u = gk + a_copy * x_var; 
        auto v = std::abs(u) - penalty * l1_regul; 
        a = 0.0;
        if (v > 0.0) {
            a = std::max(a_min,
                    std::min(a_max,
                        std::copysign(v,u) / (x_var + penalty * l2_regul)) );
        }
    }

    GLMNETPP_STRONG_INLINE
    static value_t 
    update_intercept(
            value_t& intercept,
            value_t& dlx,
            bool intr,
            value_t r_sum,
            value_t var)
    {
        auto d = 0.0; 
        if (intr) d = r_sum / var;
        if (d) { 
            intercept += d; 
            update_dlx(dlx, d, var);
        }
        return d;
    }

    GLMNETPP_STRONG_INLINE
    static bool
    check_kkt(value_t g, value_t l1_regul, value_t penalty) {
        return g > l1_regul * penalty;
    }

    /*
     * Computes strong map given a threshold.
     */
    template <class GType
            , class VPType
            , class SType
            , class SkipType>
    GLMNETPP_STRONG_INLINE
    static bool 
    compute_strong_map(
            const GType& g,
            const VPType& penalty,
            SType& strong_map,
            value_t tlam,
            SkipType skip_f) 
    {
        bool updated = false;
        for_each_with_skip(
                util::counting_iterator<index_t>(0),
                util::counting_iterator<index_t>(g.size()),
                [&](auto k) {
                    if (check_kkt(g(k), tlam, penalty(k))) {
                        strong_map[k] = true;
                        updated = true;
                    }
                }, 
                skip_f);
        return updated;
    }

    /*
     * Same as above, but allows users to perform any action upon kkt failure at index k.
     */
    template <class GType
            , class VPType
            , class SType
            , class FType
            , class SkipType>
    GLMNETPP_STRONG_INLINE
    static bool 
    compute_strong_map(
            const GType& g,
            const VPType& penalty,
            SType& strong_map,
            value_t tlam,
            FType f,
            SkipType skip_f) 
    {
        bool updated = false;
        for_each_with_skip(
                util::counting_iterator<index_t>(0),
                util::counting_iterator<index_t>(g.size()),
                [&](auto k) {
                    if (check_kkt(g(k), tlam, penalty(k))) {
                        strong_map[k] = true;
                        updated = true;
                        f(k);
                    }
                }, 
                skip_f);
        return updated;
    }

    /*
     * Computes strong map by computing the threshold based on the previous and current lambda.
     */
    template <class GType
            , class VPType
            , class SType
            , class SkipType>
    GLMNETPP_STRONG_INLINE
    static bool 
    compute_strong_map(
            const GType& g,
            const VPType& penalty,
            SType& strong_map,
            value_t beta,
            value_t lmda,
            value_t prev_lmda,
            SkipType skip_f) 
    {
        auto tlam = beta * (2.0 * lmda - prev_lmda);
        return compute_strong_map(g, penalty, strong_map, tlam, skip_f);
    }

    /*
     * Same as above, but allows users to perform any action upon kkt failure at index k.
     */
    template <class GType
            , class VPType
            , class SType
            , class FType
            , class SkipType>
    GLMNETPP_STRONG_INLINE
    static bool 
    compute_strong_map(
            const GType& g,
            const VPType& penalty,
            SType& strong_map,
            value_t beta,
            value_t lmda,
            value_t prev_lmda,
            FType f,
            SkipType skip_f) 
    {
        auto tlam = beta * (2.0 * lmda - prev_lmda);
        return compute_strong_map(g, penalty, strong_map, tlam, f, skip_f);
    }

    template <class RType, class XType>
    GLMNETPP_STRONG_INLINE
    static auto
    compute_grad(const RType& r, const XType& x) { return r.dot(x); }

public:
    GLMNETPP_STRONG_INLINE
    constexpr static bool 
    equal(value_t x, value_t y) { return x == y; }

    // TODO: may be useful to have a static interface for multi-stuff
    template <class T1, class T2>
    GLMNETPP_STRONG_INLINE
    constexpr static bool 
    equal(const Eigen::MatrixBase<T1>& x, const Eigen::MatrixBase<T2>& y) 
    {
        return (x.array() == y.array()).all(); 
    }
};

/*
 * Base class for internal implementation of any GLM elastic-net point solver.
 * This only views resources and doesn't allocate expensive data structures.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBaseViewer
    : ElnetPointInternalStaticBase<ValueType, IndexType>
{
private:
    using base_t = ElnetPointInternalStaticBase<ValueType, IndexType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using bool_t = BoolType;

    template <class IAType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalBaseViewer(
            value_t thr,
            index_t maxit,
            index_t& nin,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju
            )
        : thr_(thr)
        , maxit_(maxit)
        , nin_(nin)
        , nx_(nx)
        , mm_(nullptr, 0)
        , nlp_(nlp)
        , ia_(ia.data(), ia.size())
        , vp_(vp.data(), vp.size())
        , cl_(cl.data(), cl.rows(), cl.cols())
        , ju_(util::init_bvec<bool_t>::eval(ju))
    {}

    GLMNETPP_STRONG_INLINE void increment_passes() { ++nlp_; }
    GLMNETPP_STRONG_INLINE void coord_desc_reset() { dlx_ = 0.0; }
    GLMNETPP_STRONG_INLINE auto all_begin() const { return util::counting_iterator<index_t>(0); }
    GLMNETPP_STRONG_INLINE auto all_end() const { return util::counting_iterator<index_t>(vp_.size()); }
    GLMNETPP_STRONG_INLINE auto active_begin() const { return util::one_to_zero_iterator<index_t>(ia_.data()); }
    GLMNETPP_STRONG_INLINE auto active_end() const { return util::one_to_zero_iterator<index_t>(ia_.data() + nin_); }
    GLMNETPP_STRONG_INLINE bool is_active(index_t j) const { return mm_(j); }
    GLMNETPP_STRONG_INLINE bool has_converged() const { return dlx_ < thr_; }
    GLMNETPP_STRONG_INLINE bool has_reached_max_passes() const { return nlp_ > maxit_; }
    GLMNETPP_STRONG_INLINE bool has_reached_max_active() const { return nin_ > nx_; }
    GLMNETPP_STRONG_INLINE auto n_active() const { return nin_; }

    GLMNETPP_STRONG_INLINE
    void update_active(index_t k)
    {
        ++nin_;
        if (has_reached_max_active()) {
            throw util::max_active_reached_error();
        }
        mm_(k) = nin_; 
        ia_(nin_-1) = k+1;
    }

protected:
    using base_t::update_dlx;

    GLMNETPP_STRONG_INLINE auto& passes() { return nlp_; }
    GLMNETPP_STRONG_INLINE auto& convg_measure() { return dlx_; }
    GLMNETPP_STRONG_INLINE void set_thresh(value_t t) { thr_ = t; }
    GLMNETPP_STRONG_INLINE auto thresh() const { return thr_; }
    GLMNETPP_STRONG_INLINE const auto& endpts() const { return cl_; }
    GLMNETPP_STRONG_INLINE const auto& penalty() const { return vp_; }
    GLMNETPP_STRONG_INLINE const auto& exclusion() const { return ju_; }
    GLMNETPP_STRONG_INLINE auto active_idx(index_t k) const { return mm_(k)-1; }

    /*
     * Derived class must call this to ensure the viewers are all initialized correctly.
     */
    template <class MMType>
    GLMNETPP_STRONG_INLINE
    void construct(MMType&& mm) {
        new (&mm_) Eigen::Map<ivec_t>(mm.data(), mm.size()); 
    }

    GLMNETPP_STRONG_INLINE
    void update_dlx(value_t beta_diff, value_t x_var) {
        base_t::update_dlx(dlx_, beta_diff, x_var);
    }

    // TODO: this is only needed for multi-stuff
    template <class T>
    GLMNETPP_STRONG_INLINE 
    void update_dlx(const Eigen::MatrixBase<T>& beta_diff, value_t x_var) {
        base_t::update_dlx(dlx_, beta_diff.array().abs().maxCoeff(), x_var);
    }
    
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;
    using sp_mat_t = Eigen::SparseMatrix<value_t>;
    using ivec_t = Eigen::Matrix<index_t, Eigen::Dynamic, 1>;
    using bvec_t = util::bvec_t<bool_t>;
    using ju_t = typename std::conditional<
                    std::is_same<bool_t, bool>::value,
                    const bvec_t&,
                    Eigen::Map<const bvec_t> >::type;

private:
    value_t dlx_ = 0.0;         // measures convergence of each coord-desc
    value_t thr_;               // threshold for convergence
    const index_t maxit_;       // max number of passes
    index_t& nin_;              // number of active variables
    const index_t nx_;          // max number of active variables
    Eigen::Map<ivec_t> mm_;     // index k is j if feature k is the jth feature active
    index_t& nlp_;                      // number of passes
    Eigen::Map<ivec_t> ia_;             // active set (important that it's 1-indexed!)
    Eigen::Map<const vec_t> vp_;        // penalties on features
    Eigen::Map<const mat_t> cl_;        // limits on each feature (2 x nvars)
    ju_t ju_;                           // exclusion type
};

template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBase
    : ElnetPointInternalBaseViewer<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBaseViewer<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::ivec_t;

public:
    template <class IAType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalBase(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju
            )
        : base_t(thr, maxit, nin_, nx, nlp, ia, vp, cl, ju)
        , mm_(vp.size())
    {
        base_t::construct(mm_);
        ia.setZero();    
        mm_.setZero();
    }

private:
    index_t nin_ = 0;
    ivec_t mm_;     // index k is j if feature k is the jth feature active
};

/*
 * Base class for all derived classes for non-linear (non-gaussian) point solver.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalNonLinearBase
    : ElnetPointInternalBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::vec_t;

public:
    template <class IAType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternalNonLinearBase(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            bool intr,
            IAType& ia,
            value_t& dev0,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju
            )
        : base_t(thr, maxit, nx, nlp, ia, vp, cl, ju)
        , ga_(vp.size())
        , ixx_(vp.size(), false)
        , intr_(intr)
        , dev0_(dev0)
    {
        ga_.setZero();
    }

    GLMNETPP_STRONG_INLINE bool is_excluded(index_t k) const { return !ixx_[k]; }
    GLMNETPP_STRONG_INLINE const auto& abs_grad() const { return ga_; }
    GLMNETPP_STRONG_INLINE auto null_deviance() const { return dev0_; }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    constexpr bool check_kkt(const PointPackType&) const { return true; }

protected:
    GLMNETPP_STRONG_INLINE auto& abs_grad() { return ga_; }
    GLMNETPP_STRONG_INLINE auto& strong_map() { return ixx_; }
    GLMNETPP_STRONG_INLINE const auto& strong_map() const { return ixx_; }
    GLMNETPP_STRONG_INLINE auto has_intercept() const { return intr_; }
    GLMNETPP_STRONG_INLINE auto& null_deviance() { return dev0_; }

    template <class MaxBetaDiffSqFType>
    GLMNETPP_STRONG_INLINE
    bool has_converged_irls(
            value_t max_intr_diff_sq,
            MaxBetaDiffSqFType max_beta_diff_sq_f) const
    {
        bool ix = max_intr_diff_sq > this->thresh();
        if (!ix) {
            for (auto it = this->active_begin(); it != this->active_end(); ++it) {
                auto k = *it;
                ix = max_beta_diff_sq_f(k) > this->thresh();
                if (ix) break;
            }
        }
        return !ix;
    }

private:
    vec_t ga_;
    std::vector<bool> ixx_; 
    const bool intr_;
    value_t& dev0_;
};

} // namespace glmnetpp
