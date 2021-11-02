#pragma once
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace glmnetpp {

struct ElnetPathGaussianBase
{
private:

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
        value_t sml0;
        value_t rsqmax0;
    };

    template <class ValueType, class IntType>
    struct PointConfigPack
    {
        using value_t = ValueType;
        using int_t = IntType;

        int_t m;
        value_t ab;
        value_t dem;
        value_t alm0;
        value_t alm;
        value_t beta;
    };

protected:

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
        using value_t = ValueType;
        using int_t = IntType;

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
        const XVType& xv;
        int_t& lmu;
        AOType& ao;
        IAType& ia;
        KinType& kin;
        RSQOType& rsqo;
        ALMOType& almo;
        int_t& nlp;
        int_t& jerr;
        SetpbFType setpb_f;
        IntParamType int_param;
    };

    template <class PackType>
    auto initialize_path(PackType&& pack) const 
    {
        using pack_t = std::decay_t<PackType>;
        using value_t = typename pack_t::value_t;
        using int_t = typename pack_t::int_t;

        auto& x = pack.x;
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

    template <class IntType
            , class PackType
            , class PathConfigPackType
            , class GType>
    auto initialize_point(
            IntType m, 
            PackType&& pack, 
            PathConfigPackType&& path_pack,
            const GType& g) const
    {
        using pack_t = std::decay_t<PackType>;
        using value_t = typename pack_t::value_t;
        using int_t = typename pack_t::int_t;

        auto flmin = pack.flmin;
        auto beta = pack.beta;
        auto ni = path_pack.ni;
        const auto& ulam = pack.ulam;
        const auto& ju = pack.ju;
        const auto& vp = pack.vp;
        auto& alm = path_pack.alm;
        auto alf = path_pack.alf;
        auto omb = path_pack.omb;
        auto setpb_f = pack.setpb_f;
        auto big = pack.int_param.big;
        auto itrace = pack.int_param.itrace;

        if (itrace) setpb_f(m);
        auto alm0 = alm;
        if (flmin >= 1.0) { alm = ulam(m); }
        else if (m > 1) { alm *= alf; }
        else if (m == 0) { alm = big; } 
        else { 
            alm0 = 0.0;
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j] == 0 || (vp(j) <= 0.0)) continue; 
                // TODO: for gaussian_naive, std::abs is redundant
                // profile to see if it's worth changing
                alm0 = std::max(alm0, std::abs(g(j)) / vp(j));
            }
            alm0 /= std::max(beta, 1e-3);
            alm = alm0 * alf;
        }
        auto dem = alm * omb; 
        auto ab = alm * beta; 

        using point_config_pack_t = PointConfigPack<value_t, int_t>;
        return point_config_pack_t{m, ab, dem, alm0, alm, beta};
    }

};

} // namespace glmnetpp
