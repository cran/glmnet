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
#include <glmnetpp_bits/elnet_point/internal/binomial_base.hpp>

namespace glmnetpp {

/*
 * Dense elastic-net point solver for Binomial multi-class group method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
    util::glm_type::binomial,
    util::mode_type<util::glm_type::binomial>::multi_class_group,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalBinomialMultiClassGroupBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialMultiClassGroupBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<ValueType, IndexType, BoolType>;
    using typename base_t::state_t;
    using typename base_t::mat_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class GType
            , class XType
            , class YType
            , class WType
            , class XVType
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
            GType& g,
            value_t& dev0,
            const XType& X,
            const YType& y,
            const WType& w,
            const XVType& xv,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(intr, thr, maxit, nx, nlp, ia, g, dev0, y, w, xv, vp, cl, ju, int_param)
        , X_(X.data(), X.rows(), X.cols())
    {
        base_t::construct(
                [&](index_t ic) { base_t::initialize_resid(ic); },
                [&](index_t j, auto& grad_buff) { return compute_abs_grad(j, grad_buff); });
    }

    template <class DiffType>
    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, const DiffType& beta_diff) {
        auto& r = this->resid();
        for (int ic = 0; ic < r.cols(); ++ic) {
            gaussian_naive_t::update_resid(
                    r.col(ic), beta_diff(ic), 
                    X_.col(k).cwiseProduct(this->weight()));
        }
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType&) {
        base_t::update_beta(k, [&](index_t k, auto& buff) { compute_grad(k, buff); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(const PointPackType& pack) {
        base_t::setup_wls(pack.l1_regul(), pack.l2_regul(),
                [&](index_t ic, value_t t) { base_t::initialize_resid(ic, t); });
        auto& r = this->resid();
        for (int ic = 0; ic < r.cols(); ++ic) {
            gaussian_naive_t::update_intercept(
                    this->intercept()(ic), r.col(ic), this->convg_measure(),
                    this->has_intercept(), r.col(ic).sum(), 1., this->weight());
        }
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        return base_t::update_irls(pack.l1_regul(), 
                [&](index_t ic, auto& sc) { 
                    std::for_each(this->active_begin(), this->active_end(),
                            [&](index_t k) {
                                sc += this->beta(k)(ic) * X_.col(k); 
                            });
                },
                [&](index_t k) { base_t::initialize_resid(k); },
                [&](index_t k, auto& grad_buff) { return compute_abs_grad(k, grad_buff); });
    }

private:
    template <class DestType>
    GLMNETPP_STRONG_INLINE
    void compute_grad(index_t k, DestType&& dest) const {
        dest.noalias() = this->resid().transpose() * X_.col(k);
    }

    template <class DestType>
    GLMNETPP_STRONG_INLINE
    value_t compute_abs_grad(index_t k, DestType&& dest) const {
        compute_grad(k, dest);
        return dest.norm();
    }

    Eigen::Map<const mat_t> X_;
};

} // namespace glmnetpp

