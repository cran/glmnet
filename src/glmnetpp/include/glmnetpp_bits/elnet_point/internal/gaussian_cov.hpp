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

    void update_active(index_t k) {
        base_t::update_active(k);

        for (index_t j = 0; j < X_.cols(); ++j) {
            if (this->is_excluded(j)) continue;

            // Note: j == k case is after the is_active(j) check in legacy code.
            // Since base_t::update_active adds mm_(k),
            // it will now be considered active,
            // so we have to first check this case.
            if (j == k) { 
                c_(j, nin_-1) = xv_(j); 
                continue; 
            }
            if (this->is_active(j)) {
                c_(j, nin_-1) = c_(k, mm_(j)-1); 
                continue;
            }
            c_(j, nin_-1) = X_.col(j).dot(X_.col(k));
        }
    }

private:
    using typename base_t::mat_t;
    using base_t::xv_;
    using base_t::c_;
    using base_t::nin_;
    using base_t::mm_;

    Eigen::Map<const mat_t> X_; // data matrix
};

} // namespace glmnetpp
