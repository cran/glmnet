#pragma once
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/macros.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <Eigen/Core>

namespace glmnetpp {

/*
 * Common base class across all Cox path-solvers.
 * Cox regression uses IRLS like Poisson, but with Cox-specific
 * deviance, gradient, and Hessian computations.
 */
struct ElnetPathCoxBase
    : ElnetPathBase
{
private:
    using base_t = ElnetPathBase;

protected:
    using state_t = util::control_flow;

    /*
     * FitPack for Cox path solver.
     * Contains survival-specific data in addition to base FitPack.
     */
    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class XType
            , class SurvType      // Survival response type (contains event times, status, etc.)
            , class GType         // Linear predictor (eta)
            , class QType         // Original weights
            , class ULamType
            , class A0Type
            , class AOType
            , class IAType
            , class KinType
            , class DevType
            , class ALMType
            , class SetpbFType
            , class IntParamType>
    struct FitPack
    {
        using sub_pack_t = typename base_t::template FitPack<
            ValueType
            , JUType
            , VPType
            , CLType
            , IntType
            , XType
            , ULamType
            , AOType
            , IAType
            , KinType
            , ALMType
            , SetpbFType
            , IntParamType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;

        int_t& err_code() const { return sub_pack.err_code(); }
        int_t path_size() const { return sub_pack.path_size(); }

        sub_pack_t sub_pack;
        const SurvType& surv;     // Survival data (event times, status, start times, etc.)
        GType& g;                 // Linear predictor (eta)
        const QType& q;           // Original weights
        A0Type& a0;               // Intercepts (all zeros for Cox)
        value_t& dev0;            // Null deviance
        DevType& dev;             // Deviance at each lambda
    };

    template <class FitPackType>
    auto initialize_path(const FitPackType& pack) const
    {
        auto out = base_t::initialize_path(pack.sub_pack);
        // Original Fortran coxnet1 multiplies sml by 100 (line 3654 in glmnet5dpclean.f)
        // This makes the stopping criterion less sensitive for Cox models
        out.sml *= 100;
        return out;
    }

    template <class IntType
            , class ValueType
            , class FitPackType
            , class PathConfigPackType
            , class ElnetPointType>
    auto initialize_point(
            IntType m,
            ValueType&& lmda_curr,
            const FitPackType& pack,
            const PathConfigPackType& path_pack,
            const ElnetPointType& elnet_point) const
    {
        return base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, elnet_point.abs_grad());
    }

    template <class FitPackType
            , class PointConfigPackType
            , class PathConfigPackType
            , class ElnetPointType>
    state_t process_point_fit(
            const FitPackType& pack,
            const PathConfigPackType& path_pack,
            const PointConfigPackType& point_pack,
            const ElnetPointType& elnet_point) const
    {
        using int_t = typename std::decay_t<PointConfigPackType>::int_t;

        auto& sp = pack.sub_pack;
        auto& a = sp.ao;
        auto& a0 = pack.a0;
        auto& dev = pack.dev;
        auto mnl = path_pack.mnl;
        int_t m = point_pack.m;
        int_t nin = elnet_point.n_active();

        base_t::store_beta_compressed(
                elnet_point.active_begin(), elnet_point.active_end(),
                a.col(m), [&](int_t k) { return elnet_point.beta(k); });

        // Cox has no intercept
        a0(m) = 0.0;

        dev(m) = elnet_point.deviance();
        int_t me = (a.col(m).head(nin).array() != 0.0).count();
        auto prev_dev = (m+1 >= mnl) ? dev(m-mnl+1) : 0;
        state_t state = base_t::process_point_fit(
                m, nin, me, (dev(m)-prev_dev)/dev(m), dev(m), sp, path_pack, point_pack);
        if (state == state_t::continue_ ||
            state == state_t::break_) return state;

        return state_t::noop_;
    }

    template <class FitPackType, class ElnetPointType>
    void process_path_fit(
            const FitPackType& pack,
            const ElnetPointType& elnet_point) const
    {
        pack.g = elnet_point.prediction();
    }
};

} // namespace glmnetpp
