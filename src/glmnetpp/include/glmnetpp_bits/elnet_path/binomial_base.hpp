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
 * Common base class across all binomial path-solvers.
 */
struct ElnetPathBinomialBase
    : ElnetPathBase
{
private:
    using base_t = ElnetPathBase;

protected:
    using state_t = util::control_flow;

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class XType
            , class YType
            , class GType
            , class WType
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
        const YType& y;
        GType& g;
        const WType& w;
        bool isd;
        bool intr;
        int_t kopt;
        A0Type& a0;
        value_t& dev0;
        DevType& dev;
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
            , class GType>
    auto initialize_point(
            IntType m, 
            ValueType&& lmda_curr,
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            const GType& g) const
    {
        return base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, g);  
    }
};

/*
 * Common base class for all Binomial Two-Class path-solvers.
 */
struct ElnetPathBinomialTwoClassBase
    : ElnetPathBinomialBase
{
private:
    using base_t = ElnetPathBinomialBase;
    using base_t::initialize_point;

    template <class WType
            , class YType
            , class PType
            , class ValueType>
    GLMNETPP_STRONG_INLINE
    static auto dev2(
            const WType& w,
            const YType& y,
            const PType& p,
            ValueType pmin
            )
    {
        auto pmax = 1.0-pmin;
        auto s = 0.0;
        for (int i = 0; i < w.size(); ++i) {
            auto pi = std::min(std::max(pmin, p(i)), pmax);
            s -= w(i) * (y(i) * std::log(pi) 
                        + (1.0-y(i)) * std::log(1.0-pi));
        }
        return s;
    }

protected:
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
        return base_t::initialize_point(m, lmda_curr, pack, path_pack, elnet_point.abs_grad());  
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
        const auto& int_param = sp.int_param;

        const auto& w = pack.w;
        const auto& y = pack.y;
        auto& a0 = pack.a0;
        auto& dev = pack.dev;
        int_t m = point_pack.m;
        int_t nin = elnet_point.n_active();
        const auto& q = elnet_point.q();
        auto dev1 = elnet_point.deviance();
        auto dev0 = elnet_point.null_deviance();

        base_t::store_beta_compressed(
                elnet_point.active_begin(), elnet_point.active_end(),
                a.col(m), [&](int_t k) { return elnet_point.beta(k); } );
        a0(m) = elnet_point.intercept(); 
        auto devi = dev2(w,y,q,int_param.pmin);
        dev(m) = (dev1-devi)/dev0; 
        int_t me = (a.col(m).head(nin).array() != 0.0).array().count();
        auto prev_dev = (m == 0) ? 0 : dev(m-1); 
        state_t state = base_t::process_point_fit(
                m, nin, me, dev(m)-prev_dev, dev(m), sp, path_pack, point_pack);

        if (state == state_t::continue_ || 
            state == state_t::break_) return state;

        if (elnet_point.is_total_var_too_small()) return state_t::break_;

        return state_t::noop_;
    }

    template <class FitPackType, class ElnetPointType>
    void process_path_fit(
            const FitPackType& pack, 
            const ElnetPointType& elnet_point) const 
    {
        auto& g = pack.g;
        const auto& q = elnet_point.q();
        g.array() = (q.array()/(1.0-q.array())).log();
    }
};

/*
 * Common base class for all Binomial Multi-Class path-solvers.
 * TODO: refactor more after integrating sparse multi-class.
 */
struct ElnetPathBinomialMultiClassBase
    : ElnetPathBinomialBase
{
private:
    template <class AType
            , class MType
            , class IntType
            , class ISType>
    static auto nintot(
            const AType& a,
            const MType& m,
            IntType nin,
            ISType&& is
            ) 
    {
        auto nc = a.cols();
        is.setZero(); 
        int out = 0;
        for (int ic = 0; ic < nc; ++ic) {
            for (int j = 0; j < nin; ++j) {
                auto k = m(j)-1;
                if (is(k)) continue;
                if (!a(j,ic)) continue;
                is(k) = k+1;
                ++out;
            }
        }
        return out;
    }

protected:
    using base_t = ElnetPathBinomialBase;

    template <class ValueType, class IntType>
    struct PathConfigPack
    {
        using sub_pack_t = typename base_t::template PathConfigPack<ValueType, IntType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;
        using ivec_t = Eigen::Matrix<int_t, Eigen::Dynamic, 1>;

        sub_pack_t sub_pack;
        ivec_t is;   // common buffer for both point and path fitting
    };

    template <class ValueType, class IntType>
    struct PointConfigPack
    {
        using sub_pack_t = typename base_t::template PointConfigPack<ValueType, IntType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;
        using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;

        auto elastic_prop() const { return sub_pack.elastic_prop(); }
        auto lmda() const { return sub_pack.lmda(); }
        auto prev_lmda() const { return sub_pack.prev_lmda(); }
        auto l1_regul() const { return sub_pack.l1_regul(); }
        auto l2_regul() const { return sub_pack.l2_regul(); }

        sub_pack_t sub_pack;
        Eigen::Map<mat_t> a_slice;
    };

    template <class FitPackType>
    auto initialize_path(const FitPackType& pack) const 
    {
        using value_t = typename std::decay_t<FitPackType>::value_t;
        using int_t = typename std::decay_t<FitPackType>::int_t;
        using ivec_t = Eigen::Matrix<int_t, Eigen::Dynamic, 1>;
        auto&& sp = base_t::initialize_path(pack); 
        const auto& x = pack.sub_pack.x;
        const auto& y = pack.y;
        ivec_t is(std::max(x.cols(), y.cols()));
        using pack_config_t = PathConfigPack<value_t, int_t>;
        return pack_config_t{{std::move(sp)}, std::move(is)};
    }

    template <class IntType
            , class ValueType
            , class FitPackType
            , class PathConfigPackType
            , class GType>
    auto initialize_point(
            IntType m, 
            ValueType&& lmda_curr,
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            const GType& g) const
    {
        using pack_t = typename std::decay<FitPackType>::type;
        using value_t = typename pack_t::value_t;
        using int_t = typename pack_t::int_t;
        using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;

        auto nx = pack.sub_pack.nx;
        auto& a = pack.sub_pack.ao;
        auto nc = pack.y.cols();

        // set the new slice for coefficient storage
        Eigen::Map<mat_t> a_slice(
                a.data() + nx * nc * m, 
                nx, nc);

        auto&& sp = base_t::initialize_point(m, lmda_curr, pack, path_pack.sub_pack, g);

        using point_config_pack_t = PointConfigPack<value_t, int_t>;
        return point_config_pack_t{{std::move(sp)}, a_slice};
    }

    /*
     * Helper function that generalizes point fit processing
     * depending on how we compute the total number of nonzero entries in the active set.
     */
    template <class FitPackType
            , class PointConfigPackType
            , class PathConfigPackType
            , class ElnetPointType
            , class ActiveNznFType>
    state_t process_point_fit(
            const FitPackType& pack, 
            PathConfigPackType&& path_pack,
            PointConfigPackType&& point_pack,
            const ElnetPointType& elnet_point,
            ActiveNznFType active_nzn_f) const
    {
        using int_t = typename std::decay_t<PointConfigPackType>::int_t;

        auto& sp = pack.sub_pack;

        auto no = sp.x.rows();
        const auto& y = pack.y;
        const auto& w = pack.w;
        auto& a0 = pack.a0;
        auto& dev = pack.dev;
        int_t m = point_pack.sub_pack.m;
        auto& a_slice = point_pack.a_slice;

        auto nc = elnet_point.n_classes();
        int_t nin = elnet_point.n_active();
        const auto& q = elnet_point.q();
        const auto& sxp = elnet_point.sxp();
        auto dev1 = elnet_point.deviance();
        auto dev0 = elnet_point.null_deviance();

        auto devi = 0.0;
        for (int_t ic = 0; ic < nc; ++ic) {
            base_t::store_beta_compressed(
                    elnet_point.active_begin(), elnet_point.active_end(),
                    a_slice.col(ic), [&](int_t k) { return elnet_point.beta(k, ic); } );
            a0(ic, m) = elnet_point.intercept(ic);
            for (int_t i = 0; i < no; ++i) {
                if (y(i,ic) <= 0) continue;
                devi -= w(i) * y(i,ic) * std::log(q(i,ic) / sxp(i));
            }
        }
        int_t me = active_nzn_f();
        dev(m) = (dev1-devi)/dev0; 
        auto prev_dev = (m <= 0) ? 0 : dev(m-1);

        state_t state = base_t::process_point_fit(
                m, nin, me, dev(m) - prev_dev, dev(m), 
                sp, path_pack.sub_pack, point_pack.sub_pack
                );
        if (state == state_t::continue_ || 
            state == state_t::break_) return state;

        if (elnet_point.has_skipped_all_classes()) return state_t::break_;

        return state_t::noop_;
    }

    template <class FitPackType
            , class PointConfigPackType
            , class PathConfigPackType
            , class ElnetPointType>
    state_t process_point_fit(
            const FitPackType& pack, 
            PathConfigPackType&& path_pack,
            PointConfigPackType&& point_pack,
            const ElnetPointType& elnet_point) const
    {
        auto& a_slice = point_pack.a_slice;
        const auto& ia = pack.sub_pack.ia;
        auto nin = elnet_point.n_active();
        auto& is = path_pack.is;
        return process_point_fit(pack, path_pack, point_pack, elnet_point,
                [&]() { return nintot(a_slice, ia, nin, is); });
    }

    template <class FitPackType, class ElnetPointType>
    void process_path_fit(
            const FitPackType& pack, 
            const ElnetPointType& elnet_point) const 
    {
        auto& g = pack.g;
        auto& q = elnet_point.q();
        auto nc = pack.y.cols();

        g.array() = q.array().log();
        for (int i = 0; i < g.rows(); ++i) {
            g.row(i).array() -= g.row(i).sum()/nc;
        }
    }
};

/*
 * Common base class for all Binomial Multi-Class with group penalty path-solvers.
 */
struct ElnetPathBinomialMultiClassGroupBase
    : ElnetPathBinomialMultiClassBase
{
private:
    using base_t = ElnetPathBinomialMultiClassBase;

protected:

    /*
     * Looks stupid, but this makes the logic consistent,
     * and therefore, simplifies the other member functions.
     */
    template <class ValueType, class IntType>
    struct PathConfigPack
    {
        using sub_pack_t = typename base_t::base_t::template PathConfigPack<ValueType, IntType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;
        using ivec_t = Eigen::Matrix<int_t, Eigen::Dynamic, 1>;

        sub_pack_t sub_pack;
    };

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class XType
            , class YType
            , class GType
            , class WType
            , class ULamType
            , class XVType
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
            , YType
            , GType
            , WType
            , ULamType
            , A0Type
            , AOType
            , IAType
            , KinType
            , DevType
            , ALMType
            , SetpbFType
            , IntParamType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;

        int_t& err_code() const { return sub_pack.err_code(); }
        int_t path_size() const { return sub_pack.path_size(); }

        sub_pack_t sub_pack;
        const XVType& xv;
    };

    template <class FitPackType>
    auto initialize_path(const FitPackType& pack) const 
    {
        // Note: we want to by-pass the base class's override of PathPack
        // and go directly to its base class's version.
        auto&& sp = base_t::base_t::initialize_path(pack.sub_pack);
        using sp_t = std::decay_t<decltype(sp)>;
        using value_t = typename sp_t::value_t;
        using int_t = typename sp_t::int_t;
        using path_config_t = PathConfigPack<value_t, int_t>;
        return path_config_t{sp};
    }

    template <class IntType
            , class ValueType
            , class FitPackType
            , class PathConfigPackType
            , class GType>
    auto initialize_point(
            IntType m, 
            ValueType&& lmda_curr,
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            const GType& g) const
    {
        return base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, g);
    }

    template <class FitPackType
            , class PointConfigPackType
            , class PathConfigPackType
            , class ElnetPointType>
    state_t process_point_fit(
            const FitPackType& pack, 
            PathConfigPackType&& path_pack,
            PointConfigPackType&& point_pack,
            const ElnetPointType& elnet_point) const
    {
        auto& a_slice = point_pack.a_slice;
        auto nin = elnet_point.n_active();
        return base_t::process_point_fit(pack.sub_pack, path_pack, point_pack, elnet_point,
                [&]() { return (a_slice.col(0).head(nin).array() != 0.0).count(); });
    }

    template <class FitPackType, class ElnetPointType>
    void process_path_fit(
            const FitPackType& pack, 
            const ElnetPointType& elnet_point) const 
    {
        base_t::process_path_fit(pack.sub_pack, elnet_point);
    }
};

} // namespace glmnetpp
