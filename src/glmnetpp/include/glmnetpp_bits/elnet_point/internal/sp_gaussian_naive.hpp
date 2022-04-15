#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>
#include <Eigen/SparseCore>

namespace glmnetpp {

template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternal<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::naive,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalGaussianNaiveBase<
            ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalGaussianNaiveBase<
            ValueType, IndexType, BoolType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;

    template <class IAType
            , class YType
            , class WType
            , class XType
            , class XMType
            , class XSType
            , class XVType
            , class VPType
            , class CLType
            , class JUType>
    SpElnetPointInternal(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            YType& y,
            const WType& w,
            const XType& X,
            const XMType& xm,
            const XSType& xs,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju)
        : base_t(thr, maxit, nx, nlp, ia, xv, vp, cl, ju)
        , X_(X.rows(), X.cols(), X.nonZeros(), 
             X.outerIndexPtr(), X.innerIndexPtr(), 
             X.valuePtr(), X.innerNonZeroPtr())
        , y_(y.data(), y.size())
        , w_(w.data(), w.size())
        , xm_(xm.data(), xm.size())
        , xs_(xs.data(), xs.size())
    {
        base_t::construct([this](index_t k) { return compute_abs_grad(k); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        base_t::update_beta(k, pack.ab, pack.dem, compute_grad(k));
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        auto beta_diff_scaled = beta_diff / xs_(k);
        y_ -= beta_diff_scaled * X_.col(k);
        o_ += beta_diff_scaled * xm_(k);
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    bool check_kkt(const PointPackType& pack) {
        return base_t::check_kkt(pack.ab, [this](index_t k) { return compute_abs_grad(k); });
    }

private:
    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t k) const {
        return X_.col(k).cwiseProduct(w_).dot(
                (y_.array() + o_).matrix()
                ) / xs_(k);
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t k) const {
        return std::abs(compute_grad(k));
    }

    using typename base_t::vec_t;
    using typename base_t::mat_t;
    using spmat_t = Eigen::SparseMatrix<value_t>; 

    value_t o_ = 0.0;               // mean shift correction when updating gradient
    Eigen::Map<const spmat_t> X_;   // data matrix (sparse)
    Eigen::Map<vec_t> y_;           // unscaled residual vector
    Eigen::Map<const vec_t> w_;     // weights for each column of X_
    Eigen::Map<const vec_t> xm_;    // col-wise mean
    Eigen::Map<const vec_t> xs_;    // col-wise stddev 
                                    // (may not be actually the stddev of X, but something passed by user)
};

} // namespace glmnetpp
