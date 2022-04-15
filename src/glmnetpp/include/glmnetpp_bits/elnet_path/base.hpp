#pragma once
#include <algorithm>
#include <cmath>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/util/macros.hpp>

namespace glmnetpp {

/*
 * Common CRTP base class for all path-solvers.
 * All routines here should be those that depend on the derived class type (ElnetPathDerived),
 * i.e. those that have different implementation definitions based on the derived type.
 */
template <class ElnetPathDerived>
struct ElnetPathCRTPBase
{
private:
    using derived_t = ElnetPathDerived;
    using state_t = util::control_flow;

    // generates self()
    GLMNETPP_GENERATE_CRTP(derived_t)

public:

    /*
     * Main driver for fitting a path-wise solution to elastic net.
     * The purpose of this function is to abstract the control-flow of the fit function.
     *
     * @param   pack    the derived class FitPack.
     */
    template <class FitPackType>
    void fit(const FitPackType& pack) const
    {
        using pack_t = std::decay_t<FitPackType>;
        using value_t = typename pack_t::value_t;
        using int_t = typename pack_t::int_t;

        auto& jerr = pack.err_code();

        try {
            auto&& path_config_pack = self().initialize_path(pack);

            auto&& elnet_point = self().get_elnet_point(pack, path_config_pack);
        
            value_t lmda_curr = 0; // this makes the math work out in the point solver

            for (int_t m = 0; m < pack.path_size(); ++m) {

                auto&& point_config_pack = 
                    self().initialize_point(m, lmda_curr, pack, path_config_pack, elnet_point);
                
                try {
                    elnet_point.fit(point_config_pack);
                } 
                catch (const util::maxit_reached_error& e) {
                    jerr = e.err_code(m);
                    return;
                } 
                catch (const util::bnorm_maxit_reached_error& e) {
                    jerr = e.err_code(m);
                    return;
                }
                catch (const util::elnet_error& e) {
                    jerr = e.err_code(m);
                    break;
                } 

                state_t state = self().process_point_fit(pack, path_config_pack, point_config_pack, elnet_point);

                if (state == state_t::continue_) continue;
                if (state == state_t::break_) break;
            }

            self().process_path_fit(pack, elnet_point);
        } 
        catch (const util::elnet_error& e) {
            jerr = e.err_code(0);    
        }
    }
};

/*
 * Common base class for all path-solvers.
 * Note that these routines should NOT be inside the CRTP class,
 * as this will lead to massive code bloat.
 */
struct ElnetPathBase
{
protected:
    using state_t = util::control_flow;
    
    /*
     * Common PathConfigPack base class across every type of path-solver.
     */
    template <class ValueType, class IntType>
    struct PathConfigPack
    {
        using value_t = ValueType;
        using int_t = IntType;

        value_t omb;
        value_t alm;
        value_t alf;
        int_t ni;
        int_t mnl;
        value_t sml;
        value_t rsqmax;
    };

    /*
     * Common PointConfigPack base class across every type of path-solver.
     */
    template <class ValueType, class IntType>
    struct PointConfigPack
    {
        using value_t = ValueType;
        using int_t = IntType;

        value_t elastic_prop() const { return beta; }
        value_t lmda() const { return alm; }
        value_t prev_lmda() const { return alm0; }
        value_t l1_regul() const { return ab; }
        value_t l2_regul() const { return dem; }

        int_t m;
        value_t ab;
        value_t dem;
        value_t alm0;
        value_t alm;
        value_t beta;
    };

    /*
     * Common FitPack base class across every type of path-solver.
     */
    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class XType
            , class ULamType
            , class AOType
            , class IAType
            , class KinType
            , class ALMOType
            , class SetpbFType
            , class IntParamType>
    struct FitPack
    {
        using value_t = ValueType;
        using int_t = IntType;

        int_t& err_code() const { return jerr; }
        int_t path_size() const { return nlam; }

        value_t beta;
        const JUType& ju;
        const VPType& vp;
        const CLType& cl;
        int_t ne;
        int_t nx;
        const XType& x;
        int_t nlam;
        value_t flmin;
        const ULamType& ulam;
        value_t thr;
        int_t maxit;
        int_t& lmu;
        AOType& ao;
        IAType& ia;
        KinType& kin;
        ALMOType& almo;
        int_t& nlp;
        int_t& jerr;
        SetpbFType setpb_f;
        IntParamType int_param;
    };
                
    /*
     * Initializes common global configuration for a path-solver 
     * across every type of path-solver.
     * All outputted quantities should be those that both
     * a path-solver and a point-solver would not modify.
     */
    template <class FitPackType>
    auto initialize_path(const FitPackType& pack) const 
    {
        using pack_t = typename std::decay<FitPackType>::type;
        using value_t = typename pack_t::value_t;
        using int_t = typename pack_t::int_t;

        const auto& x = pack.x;
        auto beta = pack.beta;
        auto eps = pack.int_param.eps;
        auto mnlam = pack.int_param.mnlam;
        auto sml = pack.int_param.sml;
        auto rsqmax = pack.int_param.rsqmax;

        int_t ni = x.cols();
        value_t omb = 1.0 - beta;        

        value_t alm = 0.0; 
        value_t alf = 1.0;
        
        if (pack.flmin < 1.0) {
            auto eqs = std::max(eps, pack.flmin); 
            alf = std::pow(eqs, 1.0 / (pack.nlam - 1.));
        } 
        pack.nlp = 0;
        auto mnl = std::min(mnlam, pack.nlam);

        using pack_config_t = PathConfigPack<value_t, int_t>;
        return pack_config_t{omb, alm, alf, ni, mnl, sml, rsqmax};
    }

    /*
     * Initializes common global configuration for a point-solver
     * across every type of point-solver.
     * All outputted quantities should be those that a point-solver would not modify.
     *
     * @param   m           iteration number.
     * @param   lmda_curr   current lambda value (will be updated to next lmda after the call).
     * @param   pack        object of type FitPack of current class.
     * @param   path_pack   object of type PathConfigPack of current class.
     * @param   g           (absolute or non-absolute) gradient vector.
     */
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

        auto flmin = pack.flmin;
        auto beta = pack.beta;
        auto ni = path_pack.ni;
        const auto& ulam = pack.ulam;
        const auto& ju = pack.ju;
        const auto& vp = pack.vp;
        auto alf = path_pack.alf;
        auto omb = path_pack.omb;
        auto& setpb_f = pack.setpb_f;
        auto big = pack.int_param.big;
        auto itrace = pack.int_param.itrace;

        if (itrace) setpb_f(m);
        auto alm0 = lmda_curr;
        auto alm = alm0;
        if (flmin >= 1.0) { alm = ulam(m); }
        else if (m > 1) { alm *= alf; }
        else if (m == 0) { alm = big; } 
        else { 
            alm0 = 0.0;
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j] == 0 || (vp(j) <= 0.0)) continue; 
                alm0 = std::max(alm0, std::abs(g(j)) / vp(j));
            }
            alm0 /= std::max(beta, 1e-3);
            alm = alm0 * alf;
        }
        lmda_curr = alm;
        auto dem = alm * omb; 
        auto ab = alm * beta; 

        using point_config_pack_t = PointConfigPack<value_t, int_t>;
        return point_config_pack_t{m, ab, dem, alm0, alm, beta};
    }

    template <class IndexType, class ValueType, class FitPackType
            , class PathConfigPackType, class PointConfigPackType>
    state_t process_point_fit(
            IndexType m,
            IndexType n_active,
            IndexType me,
            ValueType prop_dev_change,
            ValueType curr_dev,
            const FitPackType& pack,
            const PathConfigPackType& path_pack,
            const PointConfigPackType& point_pack) const
    {
        auto& kin = pack.kin;
        auto& almo = pack.almo;
        auto& lmu = pack.lmu;
        auto flmin = pack.flmin;
        auto ne = pack.ne;
        auto mnl = path_pack.mnl;
        auto rsqmax = path_pack.rsqmax;
        auto sml = path_pack.sml;
        auto alm = point_pack.alm;

        kin(m) = n_active;
        almo(m) = alm; 
        lmu = m + 1;
        if (lmu < mnl || flmin >= 1.0) return state_t::continue_; 
        if ((me > ne) ||
            (prop_dev_change < sml) ||
            (curr_dev > rsqmax)) return state_t::break_; 
        return state_t::noop_;
    }

    /* Helper static routines for all derived classes */

    /*
     * Stores compressed_beta(j) = beta(k_j) where k_j = *begin at the jth iteration.
     */
    template <class Iter, class CompBetaType, class BetaType>
    GLMNETPP_STRONG_INLINE
    static void 
    store_beta_compressed(
            Iter begin,
            Iter end,
            CompBetaType&& compressed_beta,
            const BetaType& beta)
    {
        size_t j = 0;
        std::for_each(begin, end, [&](auto k) { compressed_beta(j++) = beta(k); });
    }
};

} // namespace glmnetpp
