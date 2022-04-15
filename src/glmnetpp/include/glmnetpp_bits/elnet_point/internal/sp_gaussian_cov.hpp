#pragma once
#include <glmnetpp_bits/elnet_point/internal/gaussian_cov.hpp>
#include <Eigen/SparseCore>

namespace glmnetpp {

template <class ValueType
        , class IndexType
        , class BoolType>
struct SpElnetPointInternal<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::cov,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalGaussianCovBase<
            ValueType, IndexType, BoolType>
{ 
private:
    using base_t = ElnetPointInternalGaussianCovBase<
            ValueType, IndexType, BoolType>;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;

    template <class IAType
            , class GType
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
            GType& g,
            const WType& w,
            const XType& X,
            const XMType& xm,
            const XSType& xs,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju)
        : base_t(thr, maxit, nx, nlp, ia, g, xv, vp, cl, ju)
        , w_(w.data(), w.size())
        , X_(X.rows(), X.cols(), X.nonZeros(), 
             X.outerIndexPtr(), X.innerIndexPtr(), 
             X.valuePtr(), X.innerNonZeroPtr())
        , xm_(xm.data(), xm.size())
        , xs_(xs.data(), xs.size())
    {}

    GLMNETPP_STRONG_INLINE
    void update_active(index_t k) {
        base_t::update_active(k, 
                [&](index_t j, index_t l) { 
                    return compute_sp_cov(
                            X_.col(j), X_.col(k), w_, xm_(j), xm_(k), xs_(j), xs_(k) );
                });
    }

private:
    /*
     * Computes weighted covariance between two features.
     * It is assumed that both features are weighted by w
     * and that the mean and standard deviation are weighted by w.
     */
    template <class X1Type, class X2Type, class WType>
    GLMNETPP_STRONG_INLINE
    static auto
    compute_sp_cov(
            const X1Type& x1,
            const X2Type& x2,
            const WType& w,
            value_t xm1,
            value_t xm2,
            value_t xs1,
            value_t xs2) 
    {
        auto wx2 = x2.cwiseProduct(w);
        return (x1.dot(wx2) - xm1 * xm2) / (xs1 * xs2);
    }

    using typename base_t::vec_t;
    using spmat_t = Eigen::SparseMatrix<value_t>; 

    Eigen::Map<const vec_t> w_;     // weights for each column of X_
    Eigen::Map<const spmat_t> X_;   // data matrix (sparse)
    Eigen::Map<const vec_t> xm_;    // col-wise mean
    Eigen::Map<const vec_t> xs_;    // col-wise stddev 
                                    // (may not be actually the stddev of X, but something passed by user)
};

} // namespace glmnetpp
