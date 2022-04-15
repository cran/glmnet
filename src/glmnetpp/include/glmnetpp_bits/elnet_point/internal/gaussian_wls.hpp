#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>

namespace glmnetpp {

template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
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
    using typename base_t::mat_t;

public:

    template <class XType
            , class RType
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
    ElnetPointInternal(
        value_t alm0,
        value_t almc,
        value_t alpha,
        const XType& x,
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
        : base_t(alm0, almc, alpha, r, xv, v, intr, ju, vp, cl, nx, 
                 thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp)
        , X_(x.data(), x.rows(), x.cols())
    {
        base_t::construct(
                [&](index_t j) { return compute_xv(j); },
                [&](index_t j) { return compute_abs_grad(j); });
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
    void update_resid(index_t k, value_t beta_diff) {
        gaussian_naive_t::update_resid(
                this->resid(), beta_diff,
                X_.col(k).cwiseProduct(this->weight()));
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        base_t::update_intercept(this->resid().sum());
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
        return X_.col(k).array().square().matrix().dot(this->weight());
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t k) const {
        return base_t::compute_grad(this->resid(), X_.col(k));
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t k) const {
        return std::abs(compute_grad(k));
    }

    Eigen::Map<const mat_t> X_; // data matrix
};

} // namespace glmnetpp

