#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>

namespace glmnetpp {

template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternal<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::multi,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalGaussianMultiBase<
            ValueType, IndexType, BoolType> 
{
private:
    using base_t = ElnetPointInternalGaussianMultiBase<
            ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
            ValueType, IndexType, BoolType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;

    template <class IAType
            , class YType
            , class XType
            , class WType
            , class XMType
            , class XSType
            , class XVType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    SpElnetPointInternal(
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            value_t ys0,
            YType& y,
            const XType& X,
            const WType& w,
            const XMType& xm,
            const XSType& xs,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(thr, maxit, y.cols(), nx, nlp, ia, ys0, xv, vp, cl, ju, int_param)
        , o_(y.cols())
        , X_(X.rows(), X.cols(), X.nonZeros(), 
             X.outerIndexPtr(), X.innerIndexPtr(), 
             X.valuePtr(), X.innerNonZeroPtr())
        , y_(y.data(), y.rows(), y.cols())
        , w_(w.data(), w.size())
        , xm_(xm.data(), xm.size())
        , xs_(xs.data(), xm.size())
    {
        o_.setZero();
        base_t::construct([&](index_t k, auto& g) { return compute_abs_grad(k, g); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        base_t::update_beta(k, pack.l1_regul(), pack.l2_regul(), 
                [&](index_t j, auto& g) { return compute_grad(j, g); });
    }

    template <class DiffType>
    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, const DiffType& diff) {
        for (index_t j = 0; j < y_.cols(); ++j) { 
            gaussian_naive_t::update_resid(y_.col(j), diff(j)/xs_(k), X_.col(k));
        }
        o_ += (xm_(k) / xs_(k)) * diff;
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    bool check_kkt(const PointPackType& pack) {
        return base_t::check_kkt(pack.l1_regul(), 
                [&](index_t k, auto& g) { return compute_abs_grad(k, g); });
    }

private:
    template <class GType>
    GLMNETPP_STRONG_INLINE
    void compute_grad(index_t k, GType&& g) const {
        for (index_t j = 0; j < y_.cols(); ++j) {
            // TODO: seems similar to other sparse ones...
            g(j) = X_.col(k).cwiseProduct(w_).dot(
                    (y_.col(j).array() + o_(j)).matrix()) / xs_(k);
        }
    }

    template <class GType>
    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t k, GType&& g) const {
        compute_grad(k, g);
        return g.norm();
    }

    using typename base_t::vec_t;
    using typename base_t::mat_t;
    using typename base_t::sp_mat_t;

    vec_t o_;                       // factor of w_ that is missing from the true residual 
    Eigen::Map<const sp_mat_t> X_;  // data matrix
    Eigen::Map<mat_t> y_;           // residual uncentered
    Eigen::Map<const vec_t> w_;     // weights for each datapt
    Eigen::Map<const vec_t> xm_;    // mean of each column of X 
    Eigen::Map<const vec_t> xs_;    // stddev of each column of X
};

} // namespace glmnetpp
