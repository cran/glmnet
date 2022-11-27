#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/base.hpp>
#include <glmnetpp_bits/util/types.hpp>

namespace glmnetpp {

/*
 * Base class for internal implementation of Binomial elastic-net point solver.
 * This contains all the common interface and members across all versions of binomial:
 *      - multi-class
 *      - two-class
 *      - sparse/dense
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBinomialBase
    : ElnetPointInternalNonLinearBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalNonLinearBase<ValueType, IndexType, BoolType>;

protected:
    using state_t = util::control_flow;
    using typename base_t::vec_t;
    using typename base_t::ivec_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class WType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType
            >
    ElnetPointInternalBinomialBase(
            bool isd,
            bool intr,
            index_t kopt,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            index_t no,
            index_t ni,
            value_t& dev0,
            const WType& w,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(thr, maxit, nx, nlp, intr, ia, dev0, vp, cl, ju)
        , isd_(isd)
        , kopt_(kopt)
        , pmin_(int_param.pmin)
        , vmin_((1.0+int_param.pmin) * int_param.pmin * (1.0-int_param.pmin))
        , w_(w.data(), w.size())
    {}

    constexpr GLMNETPP_STRONG_INLINE bool is_total_var_too_small() const { return xmz_ <= vmin_; }
    GLMNETPP_STRONG_INLINE auto deviance() const { return dev1_; }

protected:
    GLMNETPP_STRONG_INLINE auto optimization_type() const { return kopt_; }
    GLMNETPP_STRONG_INLINE auto& new_weight_sum() { return xmz_; }
    GLMNETPP_STRONG_INLINE auto new_weight_sum() const { return xmz_; }
    GLMNETPP_STRONG_INLINE const auto& weight() const { return w_; }
    GLMNETPP_STRONG_INLINE auto& deviance() { return dev1_; }
    GLMNETPP_STRONG_INLINE auto do_standardize() const { return isd_; }
    GLMNETPP_STRONG_INLINE auto min_prob() const { return pmin_; }

private:
    value_t xmz_ = 0.0;  // TODO: technically not needed for multi-group lasso
    const bool isd_;     // TODO: technically not needed for multi-group lasso
    const index_t kopt_; // TODO: technically not needed for multi-group lasso
    const value_t pmin_; 
    const value_t vmin_; // TODO: technically not needed for multi-group lasso
    value_t dev1_ = 0.0;
    Eigen::Map<const vec_t> w_;
};

/*
 * Base class for Binomial all uni-response methods.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBinomialUniBase
    : ElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::vec_t;

public:
    template <class IAType
            , class WType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType
            >
    ElnetPointInternalBinomialUniBase(
            bool isd,
            bool intr,
            index_t kopt,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            index_t no,
            index_t ni,
            value_t& dev0,
            const WType& w,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, no, ni, dev0, w, vp, cl, ju, int_param)
        , r_(no)
        , v_(no)
    {}

protected:
    GLMNETPP_STRONG_INLINE auto& resid() { return r_; }
    GLMNETPP_STRONG_INLINE const auto& resid() const { return r_; }
    GLMNETPP_STRONG_INLINE auto& new_weight() { return v_; }
    GLMNETPP_STRONG_INLINE const auto& new_weight() const { return v_; }

private:
    vec_t r_;
    vec_t v_;
};

/*
 * Base class for Binomial two-class methods.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBinomialTwoClassBase
    : ElnetPointInternalBinomialUniBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialUniBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::state_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class GType
            , class YType
            , class WType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternalBinomialTwoClassBase(
            bool isd,
            bool intr,
            index_t kopt,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const GType& g,
            value_t& dev0,
            const YType& y,
            const WType& w,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, y.size(), vp.size(), dev0, w, vp, cl, ju, int_param)
        , b_(vp.size() + 1)
        , xv_(vp.size())
        , bs_(vp.size() + 1)
        , q_(y.size())
        , fmax_(std::log(1.0 / int_param.pmin - 1.0))
        , fmin_(-fmax_)
        , y_(y.data(), y.size())
        , g_(g.data(), g.size())
    {}

    template <class PackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PackType& p) 
    { 
        base_t::compute_strong_map(
                this->abs_grad(), this->penalty(), this->strong_map(),
                p.elastic_prop(), p.lmda(), p.prev_lmda(), 
                [&](auto k) { return this->strong_map()[k] || !this->exclusion()[k]; });
    }

    constexpr GLMNETPP_STRONG_INLINE auto class_begin() const { return util::counting_iterator<index_t>(0); }
    constexpr GLMNETPP_STRONG_INLINE auto class_end() const { return util::counting_iterator<index_t>(1); }

    GLMNETPP_STRONG_INLINE
    void update_dlx(index_t k, value_t beta_diff) {
        base_t::update_dlx(beta_diff, xv_(k));
    }

    GLMNETPP_STRONG_INLINE value_t beta(index_t k) const { return b_(k+1); }
    GLMNETPP_STRONG_INLINE value_t intercept() const { return b_(0); }

    GLMNETPP_STRONG_INLINE const auto& q() const { return q_; }

protected:
    using typename base_t::vec_t;
    using typename base_t::ivec_t;
    using typename base_t::mat_t;

    GLMNETPP_STRONG_INLINE auto& x_var() { return xv_; }
    GLMNETPP_STRONG_INLINE const auto& offset() const { return g_; }
    GLMNETPP_STRONG_INLINE value_t& beta(index_t k) { return b_(k+1); }
    GLMNETPP_STRONG_INLINE value_t& intercept() { return b_(0); }

    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t gk, value_t l1_regul, value_t l2_regul) {
        const auto& cl = this->endpts();
        base_t::update_beta(
                beta(k), gk, xv_(k), this->penalty()(k),
                cl(0,k), cl(1,k), l1_regul, l2_regul);
    }

    GLMNETPP_STRONG_INLINE
    auto update_intercept(value_t r_sum) {
        return gaussian_naive_t::update_intercept(
                intercept(), this->resid(), this->convg_measure(), this->has_intercept(),
                r_sum, this->new_weight_sum(), this->new_weight());
    }

    GLMNETPP_STRONG_INLINE
    void setup_wls() {
        bs_(0) = b_(0); 
        std::for_each(this->active_begin(), this->active_end(),
                [&](auto k) { bs_(k+1) = b_(k+1); });
    }

    template <class PredictFType>
    GLMNETPP_STRONG_INLINE
    state_t update_irls_invariants(PredictFType predict_f) 
    {
        auto& v = this->new_weight();
        auto& r = this->resid();
        const auto& w = this->weight();

        // computes the linear prediction fi at x_i 
        // and computes the corresponding probability prediction.
        // To make sure coefficients don't wander off, we set the probability to 0 or 1
        // depending on how negative/positive fi is.
        for (index_t i = 0; i < q_.size(); ++i) {
            auto fi = predict_f(i);
            if (fi < fmin_) { q_(i) = 0.0; } 
            else if (fi > fmax_) { q_(i) = 1.0; }
            else { q_(i) = 1.0/(1.0 + std::exp(-fi)); }
        }

        v.array() = w.array() * q_.array() * (1.0 - q_.array()); 
        this->new_weight_sum() = v.sum(); 
        if (this->is_total_var_too_small()) return state_t::break_; 
        r.array() = w.array() * (y_-q_).array();
        return state_t::noop_;
    }

    template <class GradFType>
    GLMNETPP_STRONG_INLINE
    state_t update_irls_strong_set(
            GradFType grad_f,
            value_t l1_regul) 
    {
        value_t diff0 = b_(0) - bs_(0);
        bool ix = base_t::has_converged_irls(
                this->new_weight_sum() * diff0 * diff0,
                [&](index_t k) { auto d = b_(k+1)-bs_(k+1); return xv_(k)*d*d; });
        if (ix) {
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

    // TODO: can this be generalized with the multi-class version?
    template <class XVFType, class GradFType>
    void construct(XVFType xv_f, GradFType grad_f) 
    {
        b_.setZero(); 
        bs_.setZero();
        
        auto& v = this->new_weight();
        auto& r = this->resid();
        const auto& w = this->weight();
        auto& ga = this->abs_grad();
        const auto& ju = this->exclusion();
        auto& xmz = this->new_weight_sum();
        auto& dev = this->deviance();
        auto& dev0 = this->null_deviance();

        auto q0 = w.dot(y_);
        if (q0 <= this->min_prob()) {
            throw util::prob_min_reached_error(0);
        }
        if (q0 >= 1.-this->min_prob()) {
            throw util::prob_max_reached_error(0);
        }

        if (!this->has_intercept()) q0 = 0.5;
        auto bz = 0.0; 
        if (this->has_intercept()) bz = std::log(q0/(1.0-q0));

        dev = 0.0;
        xmz = 0.0;
        if ((g_.array() == 0).all()) {
            auto vi = q0*(1.-q0); 
            b_(0) = bz; 
            v = vi * w;
            r.array() = w.array() * (y_.array()-q0); 
            q_.array() = q0; 
            xmz = vi; 
            dev = -(bz*q0 + std::log(1.0-q0));
        } 
        else {
            b_(0) = 0.0;
            if (this->has_intercept()) { 
                b_(0) = azero(y_,g_,w); 
            }
            q_.array() = 1.0 / (1.0 + (-b_(0)-g_.array()).exp()); 
            v.array() = w.array() * q_.array() * (1.0-q_.array()); 
            r.array() = w.array() * (y_-q_).array(); 
            xmz = v.sum();
            dev = -(b_(0)*q0 + 
                    w.dot( (y_.array()*g_.array() + 
                            (1.0-q_.array()).log()).matrix() ));
        }

        // if we approximate the Hessian with an upper bound for approximate Newton algorithms
        if (this->optimization_type() > 0) {
            if (this->do_standardize() && this->has_intercept()) { xv_.array() = 0.25; }
            else { 
                for (index_t j = 0; j < xv_.size(); ++j) {
                    if (ju[j]) {
                        xv_(j) = 0.25 * xv_f(j); 
                    }
                }
            }
        }

        dev0 = dev;
        for (index_t i = 0; i < y_.size(); ++i) {
            if (y_(i) > 0.0) dev0 += w(i)*y_(i)*std::log(y_(i));
            if (y_(i) < 1.0) dev0 += w(i)*(1.0-y_(i))*std::log(1.0-y_(i));
        }
        this->set_thresh(this->thresh() * dev0);

        for (index_t j = 0; j < ga.size(); ++j) {
            if (!ju[j]) continue; 
            ga(j) = std::abs(grad_f(j));
        }
    }

private:
    template <class YType
            , class GType
            , class QType>
    GLMNETPP_STRONG_INLINE 
    static auto azero(
            const YType& y,
            const GType& g,
            const QType& q
            )
    {
        auto n = y.size();
        Eigen::VectorXd e(n);
        Eigen::VectorXd p(n);
        Eigen::VectorXd w(n);
        auto az = 0.0;
        e.array() = (-g).array().exp();
        auto qy = q.dot(y);
        p.array() = 1./(1. + e.array());
        while (1) {
            w.array() = q.array() * p.array() * (1.0 - p.array());
            auto d = (qy - q.dot(p)) / w.sum();
            az += d;
            if (std::abs(d) < 1e-7) break;
            auto ea0 = std::exp(-az);
            p.array() = 1./(1. + ea0 * e.array());
        }
        return az;
    }

    vec_t b_;                       // coefficients + intercept at index 0
    vec_t xv_;                      // new weighted variance of columns of x
    vec_t bs_;                      // old coefficients + intercept at index 0
    vec_t q_;                       // probability predictions
    const value_t fmax_;            // max linear prediction
    const value_t fmin_;            // min linear prediction
    Eigen::Map<const vec_t> y_;     // original y response
    Eigen::Map<const vec_t> g_;     // offsets
};

// ========================================================================
// Dense Multi-class Base classes
// ========================================================================

/*
 * Base class for Binomial all multi-response methods.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBinomialMultiBase
    : ElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::state_t;
    using typename base_t::mat_t;
    using typename base_t::vec_t;

public:
    template <class IAType
            , class GType
            , class YType
            , class WType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternalBinomialMultiBase(
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
            const YType& y,
            const WType& w,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, y.rows(), vp.size(), dev0, w, vp, cl, ju, int_param)
        , nc_(y.cols())
        , exmx_(int_param.exmx)
        , exmn_(-exmx_)
        , emin_(int_param.pmin / (1.0 - int_param.pmin))
        , emax_(1.0 / emin_)
        , b_(vp.size() + 1, y.cols())
        , bs_(vp.size() + 1, y.cols())
        , q_(y.rows(), y.cols())
        , sxp_(y.rows())
        , y_(y.data(), y.rows(), y.cols())
        , g_(g.data(), g.rows(), g.cols())
    {
        b_.setZero();
        bs_.setZero();
        sxp_.setZero();
    }

    GLMNETPP_STRONG_INLINE const auto& q() const { return q_; }
    GLMNETPP_STRONG_INLINE const auto& sxp() const { return sxp_; }
    GLMNETPP_STRONG_INLINE auto class_begin() const { return util::counting_iterator<index_t>(0); }
    GLMNETPP_STRONG_INLINE auto class_end() const { return util::counting_iterator<index_t>(nc_); }
    GLMNETPP_STRONG_INLINE bool is_excluded(index_t k) const { return !this->strong_map()[k]; }
    GLMNETPP_STRONG_INLINE auto n_classes() const { return nc_; }
    GLMNETPP_STRONG_INLINE value_t beta(index_t k, index_t ic) const { return b_(k+1, ic); }
    GLMNETPP_STRONG_INLINE value_t intercept(index_t ic) const { return b_(0, ic); }

protected:
    GLMNETPP_STRONG_INLINE auto& q() { return q_; }
    GLMNETPP_STRONG_INLINE value_t log_mean_pred_max() const { return exmx_; }
    GLMNETPP_STRONG_INLINE value_t log_mean_pred_min() const { return exmn_; }
    GLMNETPP_STRONG_INLINE value_t mean_min() const { return emin_; }
    GLMNETPP_STRONG_INLINE value_t mean_max() const { return emax_; }
    GLMNETPP_STRONG_INLINE const auto& y() const { return y_; }
    GLMNETPP_STRONG_INLINE auto& sxp() { return sxp_; }
    GLMNETPP_STRONG_INLINE const auto& offset() const { return g_; }
    GLMNETPP_STRONG_INLINE auto& beta() { return b_; }
    GLMNETPP_STRONG_INLINE const auto& beta() const { return b_; }
    GLMNETPP_STRONG_INLINE auto& old_beta() { return bs_; }
    GLMNETPP_STRONG_INLINE const auto& old_beta() const { return bs_; }

    // Helper routine to initialize residual.
    // This is specific to multi-class binomial.
    template <class RType, class YT, class VT>
    GLMNETPP_STRONG_INLINE
    void initialize_resid(
            RType&& r,
            const Eigen::MatrixBase<YT>& y,
            const Eigen::MatrixBase<VT>& v) {
        r = this->weight().cwiseProduct(y - v);
    }
    template <class RType, class YT, class VT>
    GLMNETPP_STRONG_INLINE
    void initialize_resid(
            RType&& r,
            const Eigen::MatrixBase<YT>& y,
            const Eigen::MatrixBase<VT>& v,
            value_t scale) {
        r = this->weight().cwiseProduct(y - v) / scale;
    }
    template <class RType>
    GLMNETPP_STRONG_INLINE
    void initialize_resid(
            index_t ic,
            RType&& r) {
        initialize_resid(r, y_.col(ic), q_.col(ic).cwiseQuotient(sxp_));
    }
    template <class RType>
    GLMNETPP_STRONG_INLINE
    void initialize_resid(
            index_t ic,
            RType&& r,
            value_t scale) {
        initialize_resid(r, y_.col(ic), q_.col(ic).cwiseQuotient(sxp_), scale);
    }

    template <class PackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PackType& p) 
    { 
        base_t::compute_strong_map(
                this->abs_grad(), this->penalty(), this->strong_map(),
                p.elastic_prop(), p.lmda(), p.prev_lmda(), 
                [&](auto k) { return !this->is_excluded(k) || !this->exclusion()[k]; });
    }

    // TODO: this is probably generalizable with two-class.
    GLMNETPP_STRONG_INLINE
    void construct() {
        auto no = y_.rows();
        auto nc = y_.cols();
        auto& dev0 = this->null_deviance();
        auto& dev = this->deviance();
        const auto& w = this->weight();

        dev0 = 0.0;
        for (index_t ic = 0; ic < nc; ++ic) {
            auto q0 = w.dot(y_.col(ic));
            if (q0 <= this->min_prob()) {
                throw util::prob_min_reached_error(ic);
            }
            if (q0 >= 1.0 - this->min_prob()) {
                throw util::prob_max_reached_error(ic);
            }
            if (!this->has_intercept()) {
                q0 = 1.0 / nc;
                b_(0, ic) = 0.0;
            }
            else {
                b_(0, ic) = std::log(q0);
                dev -= q0 * b_(0, ic);
            }
        }

        if (!this->has_intercept()) dev = std::log(nc);

        if ((g_.array() == 0).all()) {
            b_.row(0).array() -= b_.row(0).sum() / nc;
            for (index_t ic = 0; ic < nc; ++ic) {
                q_.col(ic).array() = std::exp(b_(0, ic));
                sxp_ += q_.col(ic);
            }
        } 
        else {
            for (index_t i = 0; i < no; ++i) {
                g_.row(i).array() -= g_.row(i).sum() / nc;
            }
            if (this->has_intercept()) kazero(b_.row(0));
            dev = 0.0;
            for (index_t ic = 0; ic < nc; ++ic) {
                q_.col(ic).array() = b_(0,ic) + g_.col(ic).array();
                dev -= w.dot( (y_.col(ic).array() * q_.col(ic).array()).matrix() );
                q_.col(ic).array() = q_.col(ic).array().exp();
                sxp_ += q_.col(ic);
            }
            vec_t sxpl = (w.array() * sxp_.array().log()).matrix();
            for (index_t ic = 0; ic < nc; ++ic) {
                dev += y_.col(ic).dot(sxpl);
            }
        }

        for (index_t ic = 0; ic < nc; ++ic) {
            for (index_t i = 0; i < no; ++i) {
                if (y_(i,ic) > 0) dev0 += w(i) * y_(i,ic) * std::log(y_(i,ic));
            }
        }
        dev0 += dev;
        this->set_thresh(this->thresh() * dev0);
    }

    GLMNETPP_STRONG_INLINE
    bool update_strong_map(value_t l1_regul)
    {
        return base_t::compute_strong_map(
            this->abs_grad(), this->penalty(), this->strong_map(), l1_regul,
            [&](auto k) { return this->strong_map()[k] || !this->exclusion()[k]; });
    }

    template <class PredBuffType
            , class OffsetType
            , class QType
            , class UpdatePredictionFType>
    GLMNETPP_STRONG_INLINE
    void update_irls_class(
            PredBuffType&& pred_buff,
            value_t intr,
            const OffsetType& offset,
            QType&& q,
            UpdatePredictionFType update_prediction_f) 
    {
        pred_buff.array() = intr + offset.array();
        update_prediction_f(pred_buff);
        pred_buff.array() = pred_buff.array().max(this->log_mean_pred_min()).min(this->log_mean_pred_max());
        this->sxp() -= q;
        q.array() = (this->mean_min() * this->sxp().array()).max(
                        pred_buff.array().exp()).min(
                            this->mean_max() * this->sxp().array());
        this->sxp() += q;
    }

private:
    template <class AZType>
    GLMNETPP_STRONG_INLINE
    auto kazero(AZType&& az)
    {
        const auto& w = this->weight();
        az.setZero();
        mat_t e = this->offset().array().exp().matrix();
        vec_t s = e.rowwise().sum();
        double dm;
        auto n = this->y().rows();
        auto kk = this->y().cols();
        do {
            dm = 0.0;
            for (index_t k = 0; k < kk; ++k) {
                auto t = 0.0;
                auto u = 0.0;
                for (index_t i = 0; i < n; ++i) {
                    auto pik = e(i,k)/s(i);
                    t += w(i) * (y_(i,k) - pik);
                    u += w(i) * pik * (1.0-pik);
                }
                auto d = t/u;
                az(k) += d;
                auto ed = std::exp(d);
                dm = std::max(dm, std::abs(d));
                for (index_t i = 0; i < n; ++i) {
                    auto z = e(i,k);
                    e(i,k) = z * ed;
                    s(i) += -z + e(i,k);
                }
            }
        } while(dm >= 1e-7);
        az.array() -= az.sum() / kk;
    }

    const index_t nc_;      // number of classes
    const value_t exmx_;    // max linear prediction
    const value_t exmn_;    // min linear prediction
    const value_t emin_;    // min probability prediction
    const value_t emax_;    // max probability prediction

    mat_t b_;               // matrix of coefficients with intercepts at row 0
    mat_t bs_;              // matrix of old coefficients with intercepts at row 0
    mat_t q_;               // matrix of probability predictions
    vec_t sxp_;             // sum of exponential terms to normalize the probabilities
    Eigen::Map<const mat_t> y_; // original y response
    Eigen::Map<mat_t> g_;       // offsets
};

/*
 * Base class for Binomial multi-class (non-group lasso) methods.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBinomialMultiClassBase
    : ElnetPointInternalBinomialMultiBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialMultiBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::state_t;
    using typename base_t::vec_t;
    using typename base_t::ivec_t;
    using typename base_t::mat_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class GType
            , class YType
            , class WType
            , class VPType
            , class CLType
            , class JUType
            , class ISType
            , class IntParamType>
    ElnetPointInternalBinomialMultiClassBase(
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
            const YType& y,
            const WType& w,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            ISType& is,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, g, dev0, y, w, vp, cl, ju, int_param)
        , xv_(vp.size(), y.cols())
        , di_(y.rows())
        , r_(y.rows())
        , v_(y.rows())
        , pfm_((1.0 + int_param.pmin) * int_param.pmin)
        , pfx_((1.0 - int_param.pmin) * (1.0 - int_param.pmin)) 
        , is_(is.data(), is.size())
        , bs_ic_(nullptr, 0)
        , b_ic_(nullptr, 0)
        , q_ic_(nullptr, 0)
        , y_ic_(nullptr, 0)
        , xv_ic_(nullptr, 0)
        , g_ic_(nullptr, 0)
    {}

    using base_t::beta;

    template <class PackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PackType& p) 
    { 
        base_t::initialize(p);
        ig_ = false;
    }

    GLMNETPP_STRONG_INLINE bool has_skipped_all_classes() const { return !ig_; }
    GLMNETPP_STRONG_INLINE void reset_converged() { ix_ = false; }

    GLMNETPP_STRONG_INLINE
    void update_dlx(index_t k, value_t beta_diff) {
        base_t::update_dlx(beta_diff, xv_ic_(k));
    }

    GLMNETPP_STRONG_INLINE value_t& beta(index_t k) { return b_ic_(k+1); }
    GLMNETPP_STRONG_INLINE value_t beta(index_t k) const { return b_ic_(k+1); }

protected:
    using base_t::initialize_resid;

    GLMNETPP_STRONG_INLINE auto& resid() { return r_; }
    GLMNETPP_STRONG_INLINE const auto& resid() const { return r_; }
    GLMNETPP_STRONG_INLINE auto& curr_x_var() { return xv_ic_; }
    GLMNETPP_STRONG_INLINE auto curr_intercept() const { return b_ic_(0); }
    GLMNETPP_STRONG_INLINE auto& curr_intercept() { return b_ic_(0); }
    GLMNETPP_STRONG_INLINE auto& curr_q() { return q_ic_; }
    GLMNETPP_STRONG_INLINE auto& curr_offset() { return g_ic_; }
    GLMNETPP_STRONG_INLINE auto& prediction_buffer() { return di_; }
    GLMNETPP_STRONG_INLINE auto& new_weight() { return v_; }
    GLMNETPP_STRONG_INLINE const auto& new_weight() const { return v_; }

    GLMNETPP_STRONG_INLINE 
    void initialize_resid(index_t ic) {
        base_t::initialize_resid(ic, this->resid());
    }

    GLMNETPP_STRONG_INLINE
    auto update_intercept(value_t r_sum) {
        return gaussian_naive_t::update_intercept(
                curr_intercept(), this->resid(), this->convg_measure(), this->has_intercept(),
                r_sum, this->new_weight_sum(), this->new_weight());
    }

    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t gk, value_t l1_regul, value_t l2_regul) {
        const auto& cl = this->endpts();
        base_t::update_beta(
                beta(k), gk, xv_ic_(k), this->penalty()(k),
                cl(0,k), cl(1,k), l1_regul, l2_regul);
    }

    GLMNETPP_STRONG_INLINE
    state_t setup_wls(size_t ic) {
        const auto& g = this->offset();
        const auto& y = this->y();
        auto& q = this->q();
        auto& b = this->beta();
        auto& bs = this->old_beta();
        const auto& sxp = this->sxp();

        // set the viewers to current class
        new (&bs_ic_) Eigen::Map<vec_t>      (bs.col(ic).data(), bs.rows());
        new (&b_ic_ ) Eigen::Map<vec_t>      (b.col(ic).data(), b.rows());
        new (&q_ic_ ) Eigen::Map<vec_t>      (q.col(ic).data(), q.rows());
        new (&y_ic_ ) Eigen::Map<const vec_t>(y.col(ic).data(), y.rows());
        new (&xv_ic_) Eigen::Map<vec_t>      (xv_.col(ic).data(), xv_.rows());
        new (&g_ic_ ) Eigen::Map<const vec_t>(g.col(ic).data(), g.rows());

        auto& v = this->new_weight();
        auto& xmz = this->new_weight_sum();
        const auto& w = this->weight();

        // do some setup
        bs_ic_(0) = b_ic_(0); 
        std::for_each(this->active_begin(), this->active_end(),
                [&](auto i) { bs_ic_(i+1) = b_ic_(i+1); });
        xmz = 0.0;
        for (index_t i = 0; i < y_ic_.size(); ++i) {
            auto pic = q_ic_(i) / sxp(i);
            if (pic < pfm_) { pic = 0.0; v(i) = 0.0; }
            else if (pic > pfx_) { pic = 1.0; v(i) = 0.0; }
            else { v(i) = w(i) * pic * (1.0 - pic); xmz += v(i); }
            this->resid()(i) = w(i) * (y_ic_(i) - pic);
        }
        if (this->is_total_var_too_small()) return state_t::continue_;
        ig_ = true;

        return state_t::noop_;
    }

    template <class UpdateYPredFType
            , class UpdatePPredFType
            , class InitResidFType
            , class ComputeGradFType>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(
            value_t elastic_prop,
            value_t l1_regul,
            UpdateYPredFType update_y_pred_f,
            UpdatePPredFType update_p_pred_f,
            InitResidFType init_resid_f,
            ComputeGradFType compute_grad_f)
    {
        auto beta = elastic_prop;
        auto ab = l1_regul;
        auto& y = this->y();
        auto& b = this->beta();
        auto& sxp = this->sxp();
        auto& q = this->q();
        auto nc = y.cols();

        auto s = -b.row(0).sum() / nc;
        b.row(0).array() += s;
        di_.array() = s;

        // Bug fix: necessary to subtract since otherwise the next loop does 1 too many iterations.
        // Original Fortran code has a memory issue here.
        // TODO: honestly, this part of the code should not be executed when max active is reached.
        // It should just throw to the top-level caller (elnet_path) and handle there.
        auto begin = this->active_begin();
        auto end = this->active_end();
        if (this->has_reached_max_active()) --end;

        std::for_each(begin, end, 
                [&](auto l) {  
                    if (this->penalty()(l) <= 0) { s = b.row(l+1).sum()/nc; }
                    else { s = elc(beta, this->endpts().col(l), b.row(l+1)); }
                    b.row(l+1).array() -= s;
                    update_y_pred_f(l, s, di_);
                });
        update_p_pred_f(di_);

        sxp.array() *= di_.array();
        std::for_each(this->class_begin(), this->class_end(), 
                [&](index_t ic) { q.col(ic).array() *= di_.array(); });
        if (this->has_reached_max_active()) throw util::max_active_reached_error();
        if (!ig_) return state_t::break_;
        if (!has_some_class_not_converged()) {
            update_abs_grad(init_resid_f, compute_grad_f);
            ix_ = base_t::update_strong_map(ab);
            if (!has_some_class_not_converged()) return state_t::break_;
        }

        return state_t::noop_;
    }

    template <class UpdatePredictionFType>
    GLMNETPP_STRONG_INLINE
    void update_irls_class(
            UpdatePredictionFType update_prediction_f) 
    {
        base_t::update_irls_class(
                di_, this->curr_intercept(), this->curr_offset(), this->curr_q(),
                update_prediction_f);
    }

    GLMNETPP_STRONG_INLINE
    void has_converged_irls_class()
    {
        // only check for convergence of current class if all previous classes seemed to have converged.
        if (!has_some_class_not_converged()) {
            auto d = bs_ic_(0) - b_ic_(0);
            ix_ = !base_t::has_converged_irls(
                    this->new_weight_sum() * d * d,
                    [&](index_t k) { 
                        auto d = b_ic_(k+1) - bs_ic_(k+1);
                        return xv_ic_(k) * d * d;
                    });
        }
    }

    template <class XVFType
            , class InitResidFType
            , class ComputeGradFType>
    void construct(
            XVFType xv_f,
            InitResidFType init_resid_f,
            ComputeGradFType compute_grad_f)
    {
        base_t::construct();
        if (this->optimization_type() > 0) {
            if (this->do_standardize() && this->has_intercept()) { xv_.array() = 0.25; }
            else { 
                base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                        [&](index_t j) { xv_.row(j).array() = 0.25 * xv_f(j); },
                        [&](index_t j) { return !this->exclusion()[j]; });
            }
        }
        update_abs_grad(init_resid_f, compute_grad_f);
    }

private:

    GLMNETPP_STRONG_INLINE bool has_some_class_not_converged() const { return ix_; }

    template <class InitResidFType
            , class ComputeGradFType>
    GLMNETPP_STRONG_INLINE
    void update_abs_grad(
            InitResidFType init_resid_f,
            ComputeGradFType compute_grad_f)
    {
        auto& ga = this->abs_grad();
        auto skip_f = [&](auto k) {
            return this->strong_map()[k] || !this->exclusion()[k];
        };
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                [&](auto k) { ga(k) = 0.; }, skip_f);

        std::for_each(this->class_begin(), this->class_end(),
                [&](auto ic) {
                    init_resid_f(ic);
                    base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                        [&](index_t k) { 
                            ga(k) = std::max(ga(k), std::abs(compute_grad_f(k)));
                        },
                        skip_f);
                });
    }

    template <class CLType
            , class AType>
    GLMNETPP_STRONG_INLINE
    auto elc(value_t parm,
             const CLType& cl,
             const AType& a)
    {
        auto n = a.size();
        auto am = a.sum()/n;
        auto out = am;
        if (parm && (n != 2)) { 
            // Note: it is VERY important that we take head(n).
            // Otherwise, is_ is reshaped to size n, which is undefined behavior for the caller.
            is_.head(n) = Eigen::VectorXi::LinSpaced(n, 0, n-1);
            std::sort(is_.data(), is_.data() + n, 
                    [&](size_t i, size_t j) { return a[i] < a[j]; });

            if (a(is_(0)) == a(is_(n-1))) {
                out = a(0); 
            } else {
                double ad = 0.0;
                if (n % 2 == 1) { ad = a(is_(n/2)); }
                else { ad = 0.5*(a(is_(n/2))+a(is_(n/2-1))); }
                if (parm == 1.0) { 
                    out = ad; 
                } else {
                    auto b1 = std::min(am, ad); 
                    auto b2 = std::max(am, ad); 
                    auto k2=1;
                    while (a(is_(k2-1)) <= b1) { ++k2; }
                    auto k1 = k2-1; 
                    while (a(is_(k2-1)) < b2) { ++k2; }
                    auto r = parm/((1.0-parm)*n); 
                    auto is = 0; 
                    auto sm = n-2*(k1-1);
                    auto s = 0.0;
                    for (int k = k1; k < k2; ++k) {
                        sm -= 2;
                        s = r * sm + am;
                        if (s > a(is_(k-1)) && s <= a(is_(k))) {
                            is = k;
                            break;
                        }
                    }
                    if (is) {
                        out = s;
                    } else {
                        auto r2 = 2.0 * r; 
                        auto s1 = a(is_(k1-1)); 
                        auto am2 = 2.0 * am;
                        auto cri = r2 * (a.array()-s1).abs().sum() + s1*(s1-am2); 
                        out = s1;
                        for (int k = k1+1; k < k2+1; ++k) {
                            s = a(is_(k-1));
                            if (s == s1) continue;
                            auto c = r2 * (a.array()-s).abs().sum() + s*(s-am2);
                            if (c < cri) { cri = c; out = s; }
                            s1 = s;
                        }
                    }
                }
            } 
        }
        out = std::max((a.array()-cl(1)).maxCoeff(),
                std::min((a.array()-cl(0)).minCoeff(), out) );
        return out;
    }

    mat_t xv_;          // matrix of weighted variance of columns of x (note: each class has different weights)
    vec_t di_;          // buffer of length(n) used for both linear predictions and probability predictions
    vec_t r_;           // residual
    vec_t v_;           // new weights
    const value_t pfm_; // min probability predictions
    const value_t pfx_; // max probability predictions
    bool ig_ = false;   // true if all classes were skipped in WLS
    bool ix_ = false;   // true if some class did not converge

    Eigen::Map<ivec_t> is_; // a buffer used only for elc

    // current views of column for given class ic (see set_class())
    Eigen::Map<vec_t> bs_ic_;
    Eigen::Map<vec_t> b_ic_;
    Eigen::Map<vec_t> q_ic_;
    Eigen::Map<const vec_t> y_ic_;
    Eigen::Map<vec_t> xv_ic_;
    Eigen::Map<const vec_t> g_ic_;
};

/*
 * Base class for Binomial multi-class group lasso methods.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalBinomialMultiClassGroupBase
    : ElnetPointInternalBinomialMultiBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialMultiBase<ValueType, IndexType, BoolType>;
    using gaussian_multi_t = ElnetPointInternalGaussianMultiBase<ValueType, IndexType, BoolType>;

protected:
    using typename base_t::state_t;
    using typename base_t::mat_t;
    using typename base_t::vec_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class GType
            , class YType
            , class WType
            , class XVType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternalBinomialMultiClassGroupBase(
            bool intr,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            GType& g,
            value_t& dev0,
            const YType& y,
            const WType& w,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(true /* isd not used */, intr, 2 /* kopt not used */, thr, maxit, nx, nlp, ia, g, dev0, y, w, vp, cl, ju, int_param)
        , bnorm_thr_(int_param.bnorm_thr)
        , bnorm_mxit_(int_param.bnorm_mxit)
        , eps_(int_param.eps)
        , xv_(xv.data(), xv.size())
        , r_(y.rows(), y.cols())
        , g_next_(y.cols())
        , del_(y.cols())
        , isc_(y.cols())
        , sc_(y.rows())
    {}

    using base_t::initialize;
    using base_t::beta;
    using base_t::intercept;

    constexpr GLMNETPP_STRONG_INLINE bool has_skipped_all_classes() const { return false; }
    GLMNETPP_STRONG_INLINE void increment_passes() { base_t::passes() += r_.cols(); }
    constexpr GLMNETPP_STRONG_INLINE void update_intercept() const {}
    GLMNETPP_STRONG_INLINE auto& beta_buffer() { return del_; }

    template <class DiffType>
    GLMNETPP_STRONG_INLINE void update_dlx(index_t k, const DiffType& diff) {
        base_t::update_dlx(diff, xv_(k));
    }

protected:
    GLMNETPP_STRONG_INLINE auto& resid() { return r_; }
    GLMNETPP_STRONG_INLINE const auto& resid() const { return r_; }
    GLMNETPP_STRONG_INLINE auto intercept() { return this->beta().row(0); }
    GLMNETPP_STRONG_INLINE auto intercept() const { return this->beta().row(0); }
    GLMNETPP_STRONG_INLINE auto beta(index_t k) { return base_t::beta().row(k+1).transpose(); }
    GLMNETPP_STRONG_INLINE auto beta(index_t k) const { return base_t::beta().row(k+1).transpose(); }

    GLMNETPP_STRONG_INLINE 
    void initialize_resid(index_t ic) {
        base_t::initialize_resid(ic, r_.col(ic));
    }
    GLMNETPP_STRONG_INLINE 
    void initialize_resid(index_t ic, value_t scale) {
        base_t::initialize_resid(ic, r_.col(ic), scale);
    }

    template <class InitResidFType
            , class ComputeAbsGradFType>
    GLMNETPP_STRONG_INLINE
    void construct(
            InitResidFType init_resid_f,
            ComputeAbsGradFType compute_abs_grad_f) 
    {
        base_t::construct();
        update_abs_grad(init_resid_f, compute_abs_grad_f);
    }

    template <class ComputeGradFType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const ComputeGradFType& compute_grad_f) {
        gaussian_multi_t::update_beta(k, beta(k), xv_(k), this->penalty()(k),
                g_next_, g_next_, l1_regul_scaled_, l2_regul_scaled_, 
                bnorm_thr_, bnorm_mxit_, isc_, 
                [&](auto i, auto) { return this->endpts()(i,k); },
                compute_grad_f);
    }

    template <class InitResidFType>
    GLMNETPP_STRONG_INLINE 
    void setup_wls(
            value_t l1_regul, 
            value_t l2_regul,
            InitResidFType init_resid_f)
    {
        const auto& q = this->q();
        const auto& sxp = this->sxp();
        value_t t = 0.0;
        for (int ic = 0; ic < q.cols(); ++ic) {
            t = std::max(t,
                    (q.col(ic).array()*
                     (1.0-q.col(ic).array()/sxp.array())/
                     sxp.array()).maxCoeff() );
        }
        if (t < eps_) throw util::below_min_variance_error();
        t *= 2.0; 
        l1_regul_scaled_ = l1_regul / t; 
        l2_regul_scaled_ = l2_regul / t;

        this->old_intercept() = this->intercept();
        std::for_each(this->active_begin(), this->active_end(),
                [&](index_t k) { this->old_beta(k).noalias() = this->beta(k); });

        for (int ic = 0; ic < r_.cols(); ++ic) {
            init_resid_f(ic, t);
        }
    }

    template <class UpdatePredictionFType
            , class InitResidFType
            , class ComputeAbsGradFType>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(
            value_t l1_regul,
            UpdatePredictionFType update_prediction_f,
            InitResidFType init_resid_f,
            ComputeAbsGradFType compute_abs_grad_f)
    {
        value_t int_diff = (this->intercept() - this->old_intercept()).array().abs().maxCoeff();
        bool ix = base_t::has_converged_irls(
                int_diff * int_diff,
                [&](index_t k) {
                    value_t b_diff = (this->beta(k) - this->old_beta(k)).array().abs().maxCoeff();
                    return b_diff * b_diff * xv_(k);
                });

        std::for_each(this->class_begin(), this->class_end(),
                [&](index_t ic) {
                    base_t::update_irls_class(
                            sc_, intercept()(ic), this->offset().col(ic), this->q().col(ic), 
                            [&](auto& sc) { update_prediction_f(ic, sc); });
                });

        intercept().array() -= intercept().sum() / this->n_classes();

        // if IRLS converged, check the strong rule if any additional features should be added.
        if (ix) {
            update_abs_grad(init_resid_f, compute_abs_grad_f);    
            ix = base_t::update_strong_map(l1_regul);
            // if no update to strong map, done.
            if (!ix) return state_t::break_;
        }

        return state_t::noop_;
    }
    
private:
    GLMNETPP_STRONG_INLINE auto old_intercept() { return base_t::old_beta().row(0); }
    GLMNETPP_STRONG_INLINE auto old_intercept() const { return base_t::old_beta().row(0); }
    GLMNETPP_STRONG_INLINE auto old_beta(index_t k) { return base_t::old_beta().row(k+1).transpose(); }
    GLMNETPP_STRONG_INLINE auto old_beta(index_t k) const { return base_t::old_beta().row(k+1).transpose(); }

    /*
     * Updates absolute gradient quantity.
     * Note that the implementation detail is slightly different from ElnetPointInternalBinomialMultiClassBase.
     * This is slightly more efficient since we can save the residual _matrix_ and not just a vector.
     */
    template <class InitResidFType
            , class ComputeAbsGradFType>
    GLMNETPP_STRONG_INLINE
    void update_abs_grad(
            InitResidFType init_resid_f,
            ComputeAbsGradFType compute_abs_grad_f)
    {
        std::for_each(this->class_begin(), this->class_end(), init_resid_f);
        auto& ga = this->abs_grad();
        auto skip_f = [&](index_t j) { return !this->is_excluded(j) || !this->exclusion()[j]; };
        // TODO: I believe this is the correct implementation
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                [&](index_t j) { ga(j) = compute_abs_grad_f(j, g_next_); },
                skip_f);
        
        // TODO: But this is how it's done in fortran...
        //base_t::for_each_with_skip(this->all_begin(), this->all_end(),
        //        [&](index_t j) { auto gj = compute_abs_grad_f(j, g_next_); ga(j) = gj*gj; },
        //        skip_f);
        //ga.array() = ga.array().sqrt();
    }

    const value_t bnorm_thr_;       // internal parameter (copied in ctor) for bnorm threshold convergence.
    const index_t bnorm_mxit_;      // internal parameter (copied in ctor) for bnorm max iteration.
    value_t eps_;                   // min variance
    value_t l1_regul_scaled_ = 0.0; // l1 regularization (lmda * elastic_prop) scaled by some variance quantity.
    value_t l2_regul_scaled_ = 0.0; // l2 regularization (lmda * elastic_prop) scaled by some variance quantity.
    Eigen::Map<const vec_t> xv_;    // weighted variance of each column of x
    mat_t r_;           // residual matrix
    vec_t g_next_;      // buffer for storing gradient
    vec_t del_;         // difference in current beta (row of matrix)
    vec_t isc_;         // buffer only needed for chkbnds
    vec_t sc_;          // buffer for storing objective values
};

// ========================================================================
// Sparse Base classes
// ========================================================================

/*
 * Base class for sparse Binomial methods.
 * This contains all extra things that sparse binomial methods require
 * (does not contain BinomialBase stuff).
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternalBinomialBase
    : ElnetPointInternalStaticBase<ValueType, IndexType>
{
private:
    using base_t = ElnetPointInternalStaticBase<ValueType, IndexType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
        ValueType, IndexType, BoolType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using bool_t = BoolType;

    template <class XType, class XBType, class XSType>
    SpElnetPointInternalBinomialBase(
                const XType& X,
                const XBType& xb,
                const XSType& xs
                )
        : X_(X.rows(), X.cols(), X.nonZeros(), 
             X.outerIndexPtr(), X.innerIndexPtr(), 
             X.valuePtr(), X.innerNonZeroPtr())
        , xb_(xb.data(), xb.size())
        , xs_(xs.data(), xs.size())
        , xm_(X.cols())
    {}

protected:

    GLMNETPP_STRONG_INLINE auto sum_weighted_resid() const { return svr_; }

    template <class NewWeightType>
    GLMNETPP_STRONG_INLINE
    void update_active(index_t k, const NewWeightType& new_weight) {
        xm_(k) = X_.col(k).dot(new_weight);
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept(value_t diff, value_t new_weight_sum) {
        if (diff) svr_ -= diff * new_weight_sum;
    }

    template <class RType, class VType>
    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, RType& r, value_t beta_diff, const VType& v, value_t new_weight_sum) {
        auto d_scaled = beta_diff/ xs_(k);
        gaussian_naive_t::update_resid(r, d_scaled, X_.col(k).cwiseProduct(v));
        o_ += d_scaled * xb_(k); 
        svr_ -= d_scaled*(xm_(k)-xb_(k)*new_weight_sum);
    }

    template <class VType>
    void update_with_new_weights(
            index_t j,
            const VType& v,
            index_t opt_type,
            value_t new_weight_sum,
            value_t& xv_j) 
    {
        xm_(j) = X_.col(j).dot(v);
        if (!opt_type) {
            xv_j = X_.col(j).cwiseProduct(X_.col(j)).dot(v);
            xv_j = (xv_j - 2.0*xb_(j)*xm_(j)+new_weight_sum*xb_(j)*xb_(j))/(xs_(j) * xs_(j));
        }
    }

    void update_shifts(value_t sum_weighted_resid) {
        svr_ = sum_weighted_resid;
        o_ = 0.0;
    }

    template <class DiType>
    void update_prediction(index_t l, value_t s, DiType& di, value_t& b0)
    {
        auto s_scaled = s/xs_(l);
        di -= s_scaled * X_.col(l);
        b0 += s_scaled * xb_(l);
    }

    template <class WType>
    GLMNETPP_STRONG_INLINE
    value_t compute_xv(index_t j, const WType& w) const {
        auto xj_sq = X_.col(j).cwiseProduct(X_.col(j));
        return (xj_sq.dot(w) - xb_(j) * xb_(j)) / (xs_(j) * xs_(j));
    }

    template <class RType, class VType>
    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t k, const RType& r, const VType& v) const {
        auto gk = X_.col(k).dot((r.array() + v.array() * o_).matrix());
        return (gk - svr_ * xb_(k))/xs_(k);
    }

private:
    using sp_mat_t = Eigen::SparseMatrix<value_t>;
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;

    value_t o_ = 0.0;    // mean shift needed in gradient update
                         // accumulates the updates coming from the centering of the columns of x.
    value_t svr_ = 0.0;  // sum of weighted residual
    Eigen::Map<const sp_mat_t> X_;  // sparse data matrix
    Eigen::Map<const vec_t> xb_;    // X column means
    Eigen::Map<const vec_t> xs_;    // X column stddevs
    vec_t xm_;  // weighted means of columns of X.
                // Note: original code includes xmz_ as the first element,
                // but in our implementation, xmz_ is a separate variable in a base class.
};

} // namespace glmnetpp
