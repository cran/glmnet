#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>

namespace glmnetpp {

template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternal<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::wls,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalGaussianWLSBase<
        ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianWLSBase<
        ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
        ValueType, IndexType, BoolType>;

protected:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::sp_mat_t;
    using typename base_t::mat_t;
    using typename base_t::vec_t;

public:

    template <class XType
            , class RType
            , class XMType
            , class XSType
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
    SpElnetPointInternal(
        value_t alm0,
        value_t almc,
        value_t alpha,
        const XType& X,
        RType& r,
        const XMType& xm,
        const XSType& xs,
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
        : base_t(alm0, almc, alpha, r, xv, v, intr, ju, vp, cl, nx, 
                 thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp)
        , X_(X.rows(), X.cols(), X.nonZeros(), 
             X.outerIndexPtr(), X.innerIndexPtr(), 
             X.valuePtr(), X.innerNonZeroPtr())
        , xm_(xm.data(), xm.size())
        , xs_(xs.data(), xs.size())
    {
        svr_ = this->resid().sum();
        base_t::construct(
                [&](index_t j) { return compute_xv(j); },
                [&](index_t j) { return std::abs(compute_abs_grad(j)); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void initialize(const PointPackType&) {
        base_t::initialize([&](index_t j) { return compute_xv(j); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType&) {
        base_t::update_beta(k, compute_grad(k));
    }

    GLMNETPP_STRONG_INLINE
    void update_rsq(index_t k, value_t diff) {
        base_t::update_rsq(k, diff);
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        auto d_scaled = beta_diff / xs_(k);
        gaussian_naive_t::update_resid(this->resid(), d_scaled,
                X_.col(k).cwiseProduct(this->weight()));
        gaussian_naive_t::update_resid(this->resid(), -d_scaled * xm_(k), this->weight());
        svr_ = this->resid().sum();
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        auto d = base_t::update_intercept(svr_);
        if (d) svr_ = this->resid().sum();
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    bool check_kkt(const PointPackType&) {
        return base_t::check_kkt(
                [&](index_t k) { return compute_xv(k); },
                [&](index_t k) { return compute_abs_grad(k); });
    }

private:
    GLMNETPP_STRONG_INLINE
    value_t compute_xv(index_t k) const {
        const auto& v = this->weight();
        value_t xv = X_.col(k).cwiseProduct(X_.col(k)).dot(v);
        xv -= 2. * xm_(k) * X_.col(k).dot(v);
        return (xv + this->new_weight_sum() * xm_(k) * xm_(k)) / (xs_(k) * xs_(k)); 
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t k) const {
        value_t d = X_.col(k).dot(this->resid());
        return (d - svr_ * xm_(k)) / xs_(k);
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t k) const {
        return std::abs(compute_grad(k));
    }

    value_t svr_ = 0.0;                 // sum of (weighted) residuals
    Eigen::Map<const sp_mat_t> X_;      // data matrix
    Eigen::Map<const vec_t> xm_;        // column means of X
    Eigen::Map<const vec_t> xs_;        // column stddevs of X
};

} // namespace glmnetpp

