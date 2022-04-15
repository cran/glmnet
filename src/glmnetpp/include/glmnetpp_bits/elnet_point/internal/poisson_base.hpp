#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/base.hpp>
#include <glmnetpp_bits/util/types.hpp>

namespace glmnetpp {

/*
 * Base class for internal implementation of Poisson elastic-net point solver.
 * This contains all the common interface and members across all versions of poisson:
 *      - naive
 *      - sparse/dense
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternalPoissonBase
    : ElnetPointInternalNonLinearBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalNonLinearBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;

protected:
    using state_t = util::control_flow;
    using typename base_t::vec_t;
    using typename base_t::ivec_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class YType
            , class GType
            , class QType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType
            >
    ElnetPointInternalPoissonBase(
            bool intr,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            index_t no,
            index_t ni,
            value_t& dev0,
            const YType& y,
            const GType& g,
            const QType& q,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(thr, maxit, nx, nlp, intr, ia, dev0, vp, cl, ju)
        , a_(ni)
        , as_(ni)
        , r_(no)
        , v_(ni)
        , w_(no)
        , fmax_(std::log(std::numeric_limits<double>::max()*0.1))
        , q_(q.data(), q.size())
        , g_(g.data(), g.size())
        , y_(y.data(), y.size())
    {
        a_.setZero();
        as_.setZero();
    }

    GLMNETPP_STRONG_INLINE auto intercept() const { return az_; }
    GLMNETPP_STRONG_INLINE auto beta(index_t k) const { return a_(k); }

    GLMNETPP_STRONG_INLINE
    void update_dlx(index_t k, value_t beta_diff) {
        base_t::update_dlx(beta_diff, v_(k));
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

protected:
    GLMNETPP_STRONG_INLINE auto& beta(index_t k) { return a_(k); }
    GLMNETPP_STRONG_INLINE auto& resid() { return r_; }
    GLMNETPP_STRONG_INLINE const auto& resid() const { return r_; }
    GLMNETPP_STRONG_INLINE auto& x_var() { return v_; }
    GLMNETPP_STRONG_INLINE const auto& x_var() const { return v_; }
    GLMNETPP_STRONG_INLINE auto new_weight_sum() const { return v0_; }
    GLMNETPP_STRONG_INLINE auto& null_deviance_intr() { return dv0_; }
    GLMNETPP_STRONG_INLINE auto null_deviance_intr() const { return dv0_; }
    GLMNETPP_STRONG_INLINE auto& intercept() { return az_; }
    GLMNETPP_STRONG_INLINE auto& weight() { return w_; }
    GLMNETPP_STRONG_INLINE const auto& weight() const { return w_; }
    GLMNETPP_STRONG_INLINE const auto& y() const { return y_; }
    GLMNETPP_STRONG_INLINE const auto& y_mean() const { return yb_; }
    GLMNETPP_STRONG_INLINE const auto& orig_weight() const { return q_; }
    GLMNETPP_STRONG_INLINE const auto& offset() const { return g_; }
    GLMNETPP_STRONG_INLINE auto max_link() const { return fmax_; }

    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, value_t gk, value_t l1_regul, value_t l2_regul) {
        const auto& cl = this->endpts();
        base_t::update_beta(
                beta(k), gk, v_(k), this->penalty()(k),
                cl(0,k), cl(1,k), l1_regul, l2_regul);
    }

    template <class UpdateFType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(UpdateFType update_f) {
        az0_ = az_; 
        std::for_each(this->active_begin(), this->active_end(),
                [&](auto k) { as_(k) = a_(k); });
        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                update_f,
                [&](auto k) { return this->is_excluded(k); });
    }

    template <class InitFType
            , class Init2FType
            , class UpdateDvrFType
            , class AbsGradFType>
    void construct(
            value_t y_mean,
            InitFType init_f,
            Init2FType init_2_f,
            UpdateDvrFType update_dvr_f,
            AbsGradFType abs_grad_f) 
    {
        bool intr = this->has_intercept();
        auto& dev0 = this->null_deviance();

        yb_ = y_mean; 
        if ((g_.array() == 0).all()) {
            if (intr) { 
                w_ = yb_ * q_;  
                az_ = std::log(yb_); 
                dv0_ = yb_ * (az_-1.0); 
                init_f(true, true);
            }
            else { 
                w_ = q_; 
                az_ = 0.;
                dv0_ = -1.0;
                init_f(true, false);
            }
        }
        else {
            w_.array() = q_.array() * (g_.array().abs().min(fmax_)).binaryExpr(g_.array(),
                        [&](auto x, auto y) { return std::copysign(x, y); }).exp();
            v0_ = w_.sum();
            if (intr) { 
                auto eaz = yb_ / v0_; 
                w_ *= eaz; 
                az_ = std::log(eaz); 
                init_f(false, true);
            }
            else { 
                az_ = 0.0; 
                init_f(false, false);
            }
        }

        v0_ = 1.0; 
        if (intr) v0_ = yb_; 
        init_2_f();
        dev0 = -yb_;
        for (int i = 0; i < y_.size(); ++i) {
            update_dvr_f(i);
        }
        dev0 -= dv0_; 
        this->set_thresh(this->thresh() * dev0); 
        gaussian_naive_t::update_abs_grad(this->abs_grad(), abs_grad_f,
                [&](index_t j) { return !this->exclusion()[j]; });
    }

    template <class GradFType>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(
            value_t l1_regul,
            GradFType grad_f) 
    {
        v0_ = w_.sum(); 
        value_t diff0 = az_ - az0_;
        bool ix = base_t::has_converged_irls(
                v0_ * diff0 * diff0,
                [&](index_t k) { auto d = a_(k)-as_(k); return v_(k)*d*d; });
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

private:
    vec_t a_;                   // coefficients
    vec_t as_;                  // old coefficients
    vec_t r_;                   // residual
    vec_t v_;                   // weighted variance of columns of x
    vec_t w_;                   // weight
    const value_t fmax_;        // max linear prediction
    value_t dv0_ = 0.0;         // null deviance
    value_t v0_ = 0.0;          // sum of weights
    value_t az_ = 0.0;          // intercept
    value_t az0_ = 0.0;         // prev intercept during WLS
    Eigen::Map<const vec_t> q_; // original weights
    Eigen::Map<const vec_t> g_; // offset
    Eigen::Map<const vec_t> y_; // y response
    value_t yb_ = 0.0;          // y mean
};

} // namespace glmnetpp
