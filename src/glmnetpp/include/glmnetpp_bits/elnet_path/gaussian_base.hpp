#pragma once
#include <algorithm>
#include <type_traits>
#include <limits>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <Eigen/Core>

namespace glmnetpp {

/*
 * Common routines across all Gaussian path-solvers.
 */
struct ElnetPathGaussianBase
    : ElnetPathBase
{
private:
    using base_t = ElnetPathBase;

protected:
    using typename base_t::state_t;
    using base_t::process_point_fit;

    /* 
     * Common FitPack base class for all Gaussian path-solvers.
     */
    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class XType
            , class ULamType
            , class XVType
            , class AOType
            , class IAType
            , class KinType
            , class RSQOType
            , class ALMOType
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
            , ALMOType
            , SetpbFType
            , IntParamType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;

        GLMNETPP_STRONG_INLINE int_t& err_code() const { return sub_pack.err_code(); }
        GLMNETPP_STRONG_INLINE int_t path_size() const { return sub_pack.path_size(); }

        sub_pack_t sub_pack;
        const XVType& xv;
        RSQOType& rsqo;
    };

    /*
     * Delegate to base class method with the base pack.
     */
    template <class FitPackType>
    GLMNETPP_STRONG_INLINE 
    auto initialize_path(const FitPackType& pack) const 
    {
        return base_t::initialize_path(pack.sub_pack);    
    }

    template <class IntType
            , class ValueType
            , class FitPackType
            , class PathConfigPackType
            , class GType>
    GLMNETPP_STRONG_INLINE 
    auto initialize_point(
            IntType m, 
            ValueType&& lmda_curr,
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            const GType& g) const
    {
        return base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, g);
    }

    /*
     * Common routine for all Gaussian path-solvers after point-solver fit.
     * See fit() in base.hpp for usage.
     *
     * @param   pack        object of FitPack of current class.
     * @param   path_pack   object of PathConfigPack of current class.
     * @param   point_pack  object of PointConfigPack of current class.
     * @param   elnet_point point-solver object.
     */
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
        using value_t = typename std::decay_t<PointConfigPackType>::value_t;

        auto& sp = pack.sub_pack;

        auto& ao = sp.ao;
        auto& rsqo = pack.rsqo;
        auto m = point_pack.m;
        auto n_active = elnet_point.n_active();
        auto rsq = elnet_point.rsq();
        auto rsq0 = elnet_point.rsq_prev();

        base_t::store_beta_compressed(
                elnet_point.active_begin(), elnet_point.active_end(),
                ao.col(m), [&](int_t k) { return elnet_point.beta(k); } );
        rsqo(m) = rsq;

        int_t me = (ao.col(m).head(n_active).array() != 0).count();
        auto prop_dev_change = (rsq == 0) ? 
            std::numeric_limits<value_t>::infinity() : 
            (rsq - rsq0) / rsq;

        state_t state = base_t::process_point_fit(
                m, n_active, me, prop_dev_change, rsq, sp, path_pack, point_pack
                );
        if (state == state_t::continue_ || 
            state == state_t::break_) return state;

        return state_t::noop_;
    }

    /*
     * Common finishing routine for all Gaussian path-solvers after (path) fit.
     * See fit() in base.hpp for usage.
     */
    template <class FitPackType, class ElnetPointType>
    GLMNETPP_STRONG_INLINE
    constexpr void process_path_fit(const FitPackType&, const ElnetPointType&) const {}
};

/*
 * Common routines across all Gaussian multi-response path-solvers.
 */
struct ElnetPathGaussianMultiBase
    : ElnetPathGaussianBase
{
private:
    using base_t = ElnetPathGaussianBase;

    template <class ValueType, class IntType>
    struct PointConfigPack
    {
        using sub_pack_t = typename base_t::template PointConfigPack<ValueType, IntType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;
        using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;

        GLMNETPP_STRONG_INLINE auto elastic_prop() const { return sub_pack.elastic_prop(); }
        GLMNETPP_STRONG_INLINE auto lmda() const { return sub_pack.lmda(); }
        GLMNETPP_STRONG_INLINE auto prev_lmda() const { return sub_pack.prev_lmda(); }
        GLMNETPP_STRONG_INLINE auto l1_regul() const { return sub_pack.l1_regul(); }
        GLMNETPP_STRONG_INLINE auto l2_regul() const { return sub_pack.l2_regul(); }

        sub_pack_t sub_pack;
        Eigen::Map<mat_t> a_slice;
    };

protected:
    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class YType
            , class IntType
            , class XType
            , class ULamType
            , class XVType
            , class AOType
            , class IAType
            , class KinType
            , class RSQOType
            , class ALMOType
            , class SetpbFType
            , class IntParamType>
    struct FitPack
    {
        using sub_pack_t = typename base_t::template FitPack<
            ValueType, JUType, VPType, CLType, IntType,
            XType, ULamType, XVType, AOType, IAType, KinType,
            RSQOType, ALMOType, SetpbFType, IntParamType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;

        int_t& err_code() const { return sub_pack.err_code(); }
        int_t path_size() const { return sub_pack.path_size(); }

        sub_pack_t sub_pack;
        YType& y;
        ValueType ys0;
    };

    template <class FitPackType>
    auto initialize_path(const FitPackType& pack) const 
    {
        return base_t::initialize_path(pack.sub_pack);
    }

    template <class IntType
            , class ValueType
            , class FitPackType
            , class PathConfigPackType
            , class ElnetPointType>
    GLMNETPP_STRONG_INLINE 
    auto initialize_point(
            IntType m, 
            ValueType&& lmda_curr,
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            const ElnetPointType& elnet_point) const
    {
        using pack_t = typename std::decay<FitPackType>::type;
        using value_t = typename pack_t::value_t;
        using int_t = typename pack_t::int_t;
        using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;

        auto nx = pack.sub_pack.sub_pack.nx;
        auto nr = pack.y.cols();
        auto& a = pack.sub_pack.sub_pack.ao;

        // set the new slice for coefficient storage
        Eigen::Map<mat_t> a_slice(
                a.data() + nx * nr * m, 
                nx, nr);

        auto&& sp = base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, elnet_point.abs_grad());

        using point_config_pack_t = PointConfigPack<value_t, int_t>;
        return point_config_pack_t{{std::move(sp)}, a_slice};
    }


    template <class FitPackType
            , class PointConfigPackType
            , class PathConfigPackType
            , class ElnetPointType>
    state_t process_point_fit(
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            PointConfigPackType&& point_pack,
            const ElnetPointType& elnet_point) const
    {
        using fit_pack_t = std::decay_t<FitPackType>;
        using value_t = typename fit_pack_t::value_t;
        using int_t = typename fit_pack_t::int_t;

        auto& point_sp = point_pack.sub_pack;
        auto nr = pack.y.cols();
        auto rsq = elnet_point.rsq();
        auto rsq0 = elnet_point.rsq_prev();
        auto n_active = elnet_point.n_active();
        auto ys0 = pack.ys0;
        auto& ao_slice = point_pack.a_slice;
        int_t m = point_pack.sub_pack.m;
        auto& rsqo = pack.sub_pack.rsqo;

        for (int j = 0; j < nr; ++j) {
            auto ao_slice_j = ao_slice.col(j);
            base_t::store_beta_compressed(
                    elnet_point.active_begin(), elnet_point.active_end(),
                    ao_slice_j, [&](auto k) { return elnet_point.beta(j, k); } );
        }
        rsqo(m) = 1.0-rsq/ys0; 
        int_t me = (ao_slice.col(0).head(n_active).array() != 0).count();
        auto prop_dev_change = (rsq == 0) ? 
            std::numeric_limits<value_t>::infinity() : 
            (rsq0 - rsq) / rsq;
        state_t state = base_t::process_point_fit(
                m, n_active, me, prop_dev_change, rsqo(m), pack.sub_pack.sub_pack, path_pack, point_sp
                );
        if (state == state_t::continue_ || 
            state == state_t::break_) return state;

        return state_t::noop_;
    }
};

} // namespace glmnetpp
