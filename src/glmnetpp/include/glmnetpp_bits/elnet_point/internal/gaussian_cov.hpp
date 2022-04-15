#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>

namespace glmnetpp {

/*
 * Dense elastic-net point solver for Gaussian covariance method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
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
            , class XType
            , class XVType
            , class VPType
            , class CLType
            , class JUType>
    ElnetPointInternal(value_t thr,
                       index_t maxit,
                       index_t nx,
                       index_t& nlp,
                       IAType& ia,
                       GType& g,
                       const XType& X,
                       const XVType& xv,
                       const VPType& vp,
                       const CLType& cl,
                       const JUType& ju)
        : base_t(thr, maxit, nx, nlp, ia, g, xv, vp, cl, ju)
        , X_(X.data(), X.rows(), X.cols())
    {}

    GLMNETPP_STRONG_INLINE
    void update_active(index_t k) {
        base_t::update_active(k, 
                [&](index_t j, index_t l) { return X_.col(j).dot(X_.col(l)); });
    }

private:
    using typename base_t::mat_t;

    Eigen::Map<const mat_t> X_; // data matrix
};

} // namespace glmnetpp
