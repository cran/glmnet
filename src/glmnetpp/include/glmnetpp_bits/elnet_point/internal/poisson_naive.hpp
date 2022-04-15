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
 * Dense elastic-net point solver for Poisson naive method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
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
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternal(
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
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(intr, thr, maxit, nx, nlp, ia, X.rows(), X.cols(), dev0, y, g, q, vp, cl, ju, int_param)
        , X_(X.data(), X.rows(), X.cols())
        , t_(X.rows())
        , f_(X.rows())
    {
        t_.array() = this->orig_weight().array() * this->y().array(); 
        this->construct(t_.sum(),
                [&](bool offset_all_zero, bool intr) {
                    if (!offset_all_zero) {
                        if (intr) {
                            this->null_deviance_intr() = 
                                t_.dot(this->offset()) - this->y_mean() * (1.0 - this->intercept());
                        } else {
                            this->null_deviance_intr() = t_.dot(this->offset()) - this->new_weight_sum(); 
                        }
                    }
                },
                [&]() { this->resid() = t_ - this->weight(); },
                [&](index_t i) { 
                    auto& dvr = this->null_deviance();
                    if (t_(i) > 0.0) dvr += t_(i) * std::log(this->y()(i));
                },
                [&](index_t j) { return compute_abs_grad(j); });
        f_.array() = this->intercept() + this->offset().array();
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        auto gk = compute_grad(k);
        base_t::update_beta(k, gk, pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        gaussian_naive_t::update_resid(
                this->resid(), beta_diff, 
                (this->weight().array() * X_.col(k).array()).matrix());
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        auto d = base_t::update_intercept(
            this->intercept(),
            this->convg_measure(),
            this->has_intercept(),
            this->resid().sum(),
            this->new_weight_sum());
        if (d) {
            gaussian_naive_t::update_resid(this->resid(), d, this->weight());
        }
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(const PointPackType&) {
        auto& v = this->x_var();
        const auto& w = this->weight();
        base_t::setup_wls([&](auto k) { v(k) = compute_xv(X_.col(k), w); });
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        auto& w = this->weight();
        const auto& q = this->orig_weight();
        auto& r = this->resid();

        f_.array() = this->intercept() + this->offset().array();
        std::for_each(this->active_begin(), this->active_end(),
                [&](index_t k) { f_ += this->beta(k) * X_.col(k); });
        w.array() = q.array() * 
            (f_.array().abs().min(this->max_link())).binaryExpr(f_.array(),
                [](auto x, auto y) { return std::copysign(x,y); }).exp();
        r=t_-w;
        return base_t::update_irls(pack.l1_regul(),
                [&](index_t k) { return compute_grad(k); });
    }

    GLMNETPP_STRONG_INLINE value_t deviance() const {
        return (t_.dot(f_) - this->new_weight_sum() - this->null_deviance_intr()) / this->null_deviance();
    }

    GLMNETPP_STRONG_INLINE const auto& prediction() const { return f_; }

private:
    using typename base_t::mat_t;
    using typename base_t::vec_t;

    // TODO: put in base class?
    template <class XType, class WType>
    GLMNETPP_STRONG_INLINE
    static value_t 
    compute_xv(const XType& x, const WType& w) { 
        return w.dot(x.array().square().matrix()); 
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t j) const {
        return base_t::compute_grad(this->resid(), X_.col(j));
    }
    
    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t j) const {
        return std::abs(base_t::compute_grad(this->resid(), X_.col(j)));
    }

    Eigen::Map<const mat_t> X_;
    vec_t t_;   // different quantity from sparse version
    vec_t f_;   // linear prediction
};

} // namespace glmnetpp
