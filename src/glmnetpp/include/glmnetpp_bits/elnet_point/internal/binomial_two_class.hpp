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
 * Dense elastic-net point solver for Binomial two-class method.
 */
template <class ValueType
        , class IndexType
        , class BoolType>
struct ElnetPointInternal<
    util::glm_type::binomial,
    util::mode_type<util::glm_type::binomial>::two_class,
    ValueType,
    IndexType,
    BoolType>
        : ElnetPointInternalBinomialTwoClassBase<ValueType, IndexType, BoolType>
{
private:
    using base_t = ElnetPointInternalBinomialTwoClassBase<ValueType, IndexType, BoolType>;
    using gaussian_naive_t = ElnetPointInternalGaussianNaiveBase<
        ValueType, IndexType, BoolType>;
    using typename base_t::state_t;

public:
    using typename base_t::value_t;
    using typename base_t::index_t;
    using typename base_t::bool_t;

    template <class IAType
            , class GType
            , class XType
            , class YType
            , class WType
            , class VPType
            , class CLType
            , class JUType
            , class IntParamType>
    ElnetPointInternal(
            bool isd,
            bool intr,
            index_t kopt,
            value_t thr,
            index_t maxit,
            index_t nx,
            index_t& nlp,
            IAType& ia,
            const GType& g,
            value_t& dev0,
            const XType& X,
            const YType& y,
            const WType& w,
            const VPType& vp,
            const CLType& cl,
            const JUType& ju,
            const IntParamType& int_param)
        : base_t(isd, intr, kopt, thr, maxit, nx, nlp, ia, g, dev0, y, w, vp, cl, ju, int_param)
        , X_(X.data(), X.rows(), X.cols())
    {
        this->construct(
                 [&](index_t j) { return compute_xv(X_.col(j), this->weight()); },
                 [&](index_t j) { return this->compute_grad(j); });
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void update_beta(index_t k, const PointPackType& pack) {
        auto gk = compute_grad(k);
        base_t::update_beta(k, gk, pack.l1_regul(), pack.l2_regul());
    }

    GLMNETPP_STRONG_INLINE
    void update_intercept() {
        base_t::update_intercept(this->resid().sum());
    }

    GLMNETPP_STRONG_INLINE
    void update_resid(index_t k, value_t beta_diff) {
        gaussian_naive_t::update_resid(
                this->resid(), beta_diff, 
                (this->new_weight().array() * X_.col(k).array()).matrix());
    }

    template <class PointPackType>
    GLMNETPP_STRONG_INLINE
    void setup_wls(const PointPackType&) {
        base_t::setup_wls();
        const auto& ixx = this->strong_map();
        auto& xv = this->x_var();
        const auto& v = this->new_weight();
        if (!this->optimization_type()) {
            base_t::for_each_with_skip(this->all_begin(), this->all_end(),
                    [&](index_t j) {
                        xv(j) = compute_xv(X_.col(j), v);
                    },
                    [&](index_t j) { return !ixx[j]; });
        }
    }

    template <class PointConfigPack>
    GLMNETPP_STRONG_INLINE
    state_t update_irls(const PointConfigPack& pack) 
    {
        auto predict_f = [&](index_t i) { 
            auto fi = this->intercept() + this->offset()(i);
            std::for_each(this->active_begin(), this->active_end(),
                    [&](auto k) { 
                        fi += this->beta(k) * X_(i,k); 
                    });
            return fi;
        };

        state_t state = base_t::update_irls_invariants(predict_f);
        if (state == state_t::break_) return state_t::break_;

        auto grad_f = [&](index_t k) { return this->compute_grad(k); };
        return base_t::update_irls_strong_set(grad_f, pack.l1_regul());
    }

private:
    using typename base_t::mat_t;

    template <class XType, class WType>
    GLMNETPP_STRONG_INLINE
    static value_t 
    compute_xv(const XType& x, const WType& w) { 
        return w.dot(x.array().square().matrix()); 
    }

    GLMNETPP_STRONG_INLINE
    value_t compute_grad(index_t j) const {
        return this->resid().dot(X_.col(j));
    }

    Eigen::Map<const mat_t> X_;
};

} // namespace glmnetpp
