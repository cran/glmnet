#pragma once
#include <cmath>
#include <utility>
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/iterator/counting_iterator.hpp>
#include <glmnetpp_bits/util/iterator/one_to_zero_iterator.hpp>
#include <glmnetpp_bits/elnet_point/internal/poisson_base.hpp>

namespace glmnetpp {

/*
 * Sparse elastic-net point solver for Poisson naive method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternal<
    util::glm_type::poisson,
    util::mode_type<util::glm_type::poisson>::naive,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalPoissonBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalPoissonBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
        ValueType, IndexType, BoolType>;
    using typename base_t::state_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class XType
            , class YType
            , class GType
            , class QType
            , class XBType
            , class XSType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    SpElnetPointInternal(
            bool intr,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            value_t& dev0,
            const XType& X,
            const YType& y,
            const GType& g,
            const QType& q,
            const XBType& xb,
            const XSType& xs,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(intr, thr, maxit, nx, nlp, ia, X.rows(), X.cols(), dev0, y, g, q, vp, cl, ju, int_param)
        , X_(X.rows(), X.cols(), X.nonZeros(), 
             X.outerIndexPtr(), X.innerIndexPtr(), 
             X.valuePtr(), X.innerNonZeroPtr())
        , xb_(xb.data(), xb.size())
        , xs_(xs.data(), xs.size())
        , t_(X.rows())
        , xm_(X.cols())
        , qy_(X.rows())
    {
        t_ = this->offset();
        qy_.array() = this->orig_weight().array() * this->y().array(); 
        this->construct(qy_.sum(),
                [&](bool offset_all_zero, bool intr) {
                    if (offset_all_zero) {
                        if (intr) {
                            uu_ = this->intercept();
                            xm_ = this->y_mean() * xb_;
                        } else {
                            xm_.setZero();
                            uu_ = 0;
                        }
                    } else {
                        if (intr) {
                            uu_ = this->intercept();
                            this->null_deviance_intr() = 
                                qy_.dot(this->offset())-this->y_mean()*(1.0-this->intercept());
                        } else {
                            uu_ = 0.0;
                            this->null_deviance_intr() = qy_.dot(this->offset()) - this->new_weight_sum(); 
                        }
                        base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                                [&](auto k) { xm_(k) = X_.col(k).dot(this->weight()); },
                                [&](auto k) { return !this->exclusion()[k]; });
                    }
                },
                [&]() {
                    auto yb = this->y_mean();
                    tt_ = yb - this->new_weight_sum() * (1.0-uu_);
                    this->resid() = qy_ - this->weight() * (1.0-uu_);
                },
                [&](index_t i) { 
                    auto& dvr = this->null_deviance();
                    if (qy_(i) > 0.0) dvr += qy_(i)*std::log(this->y()(i)); 
                },
                [&](index_t j) { return compute_abs_grad(j); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        auto gk = compute_grad(k);
        base_t::update_beta(k, gk, pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        auto d = base_t::update_intercept(
            this->intercept(),
            this->convg_measure(),
            this->has_intercept(),
            tt_ - uu_ * this->new_weight_sum(),
            this->new_weight_sum());
        uu_ += d;
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        auto d_scaled = beta_diff / xs_(k);
        gaussian_naive_t::update_resid(
                this->resid(), d_scaled, 
                X_.col(k).cwiseProduct(this->weight()));
        uu_ -= d_scaled * xb_(k);
        tt_ -= d_scaled * xm_(k);
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(const PointPackType&) {
        auto& v = this->x_var();
        const auto& w = this->weight();
        auto v0 = this->new_weight_sum();
        base_t::setup_wls(
                [&](auto k) { 
                    xm_(k) = X_.col(k).dot(w);
                    v(k) = compute_xv(X_.col(k), w, v0, xb_(k), xm_(k), xs_(k)); 
                });
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        auto& w = this->weight();
        const auto& q = this->orig_weight();
        auto& r = this->resid();

        t_ = this->offset();
        std::for_each(this->active_begin(), this->active_end(),
                [&](index_t k) { t_ += (this->beta(k) / xs_(k)) * X_.col(k); });
        w.array() = q.array() *
            ((t_.array()+uu_).abs().min(this->max_link())).matrix().binaryExpr(t_,
                [&](auto x, auto y) { return std::copysign(x,y+uu_); }).array().exp();
        r = qy_ - w*(1.0-uu_);
        tt_ = r.sum(); 

        return base_t::update_irls(pack.l1_regul(),
                [&](index_t k) {
                    xm_(k) = X_.col(k).dot(this->weight());
                    return compute_grad(k);
                });
    }

    GLMNETPP_STRONG_INLINE value_t deviance() const {
        return (qy_.dot(t_) + this->y_mean()*uu_ -this->new_weight_sum()-this->null_deviance_intr())/this->null_deviance();
    }
    GLMNETPP_STRONG_INLINE auto prediction() const { return (t_.array() + uu_).matrix(); }

private:
    using typename base_t::vec_t;
    using typename base_t::sp_mat_t;

    template <class XType, class WType>
    GLMNETPP_STRONG_INLINE
    static value_t 
    compute_xv(const XType& x,
               const WType& w,
               value_t v0,
               value_t xb,
               value_t xm,
               value_t xs) 
    {
        return (x.cwiseProduct(x).dot(w) 
                + (-2.0*xm + v0*xb)*xb ) / (xs*xs);
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t j) const {
        const auto& r = this->resid();
        auto v0 = this->new_weight_sum();
        return (X_.col(j).dot(r) - uu_*(xm_(j)-v0*xb_(j)) - xb_(j)*tt_) / xs_(j);
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t j) const {
        return std::abs(compute_grad(j));
    }
    
    value_t tt_ = 0.0;
    value_t uu_ = 0.0;
    Eigen::Map<const sp_mat_t> X_;
    Eigen::Map<const vec_t> xb_;
    Eigen::Map<const vec_t> xs_;
    vec_t t_;
    vec_t xm_;
    vec_t qy_;
};

} // namespace glmnetpp
