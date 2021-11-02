#pragma once
#include <glmnetpp_bits/elnet_point/internal/decl.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>

namespace glmnetpp {

template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
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
                       YType& y,
                       const XType& X,
                       const XVType& xv,
                       const VPType& vp,
                       const CLType& cl,
                       const JUType& ju)
        : base_t(thr, maxit, nx, nlp, ia, xv, vp, cl, ju)
        , X_(X.data(), X.rows(), X.cols())
        , y_(y.data(), y.size())
    {
        base_t::compute_abs_grad([this](index_t k) { return compute_grad(k); });
    }

    void update_beta(index_t k, value_t ab, value_t dem) {
        base_t::update_beta(k, ab, dem, compute_grad(k));
    }

    void update_resid(index_t k, value_t beta_diff) {
        y_ -= beta_diff * X_.col(k);
    }

    bool check_kkt(value_t ab) {
        return base_t::check_kkt(ab, [this](index_t k) { return compute_grad(k); });
    }

private:
    value_t compute_grad(index_t k) const {
        return y_.dot(X_.col(k));
    }

    using typename base_t::vec_t;
    using typename base_t::mat_t;

    Eigen::Map<const mat_t> X_; // data matrix
    Eigen::Map<vec_t> y_;       // scaled residual vector
                                // Note: this is slightly different from sparse version residual vector.
                                // Sparse one will not be scaled by sqrt(weights), but this one will.
};

} // namespace glmnetpp
