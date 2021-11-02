#pragma once
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_base.hpp>

namespace glmnetpp {
    
template <class SpElnetPointPolicy>
struct SpElnetPath<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::cov,
    SpElnetPointPolicy>
        : private ElnetPathGaussianBase
        , ElnetPathBase<
            SpElnetPath<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::cov,
                SpElnetPointPolicy> >
{
private:
    using gaussian_base_t = ElnetPathGaussianBase;
    using base_t = ElnetPathBase<
        SpElnetPath<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::cov,
            SpElnetPointPolicy> >;
    using elnet_point_t = SpElnetPointPolicy;

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class GType
            , class WType
            , class IntType
            , class XType
            , class ULamType
            , class XMType
            , class XSType
            , class XVType
            , class AOType
            , class IAType
            , class KinType
            , class RSQOType
            , class ALMOType
            , class SetpbFType
            , class IntParamType>
    struct FitPack :
        gaussian_base_t::template FitPack<
            ValueType, JUType, VPType, CLType, IntType,
            XType, ULamType, XVType, AOType, IAType, KinType,
            RSQOType, ALMOType, SetpbFType, IntParamType>
    {
        using base_t = gaussian_base_t::template FitPack<
            ValueType, JUType, VPType, CLType, IntType,
            XType, ULamType, XVType, AOType, IAType, KinType,
            RSQOType, ALMOType, SetpbFType, IntParamType>;
        GType& g;
        const WType& w;
        const XMType& xm;
        const XSType& xs;

        FitPack(
            ValueType beta,
            const JUType& ju,
            const VPType& vp,
            const CLType& cl,
            IntType ne,
            IntType nx,
            const XType& x,
            IntType nlam,
            ValueType flmin,
            const ULamType& ulam,
            ValueType thr,
            IntType maxit,
            const XVType& xv,
            IntType& lmu,
            AOType& ao,
            IAType& ia,
            KinType& kin,
            RSQOType& rsqo,
            ALMOType& almo,
            IntType& nlp,
            IntType& jerr,
            SetpbFType setpb_f,
            IntParamType int_param,
            GType& _g,
            const WType& _w,
            const XMType& _xm,
            const XSType& _xs
               )
            : base_t{beta, ju, vp, cl, ne, nx, x, nlam, flmin, ulam, thr, maxit, xv, 
                     lmu, ao, ia, kin, rsqo, almo, nlp, jerr, setpb_f, int_param}
            , g(_g)
            , w(_w)
            , xm(_xm)
            , xs(_xs)
        {}
    };

public:
    using gaussian_base_t::initialize_path;

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class GType
            , class WType
            , class IntType
            , class XType
            , class ULamType
            , class XMType
            , class XSType
            , class XVType
            , class AOType
            , class IAType
            , class KinType
            , class RSQOType
            , class ALMOType
            , class SetpbFType
            , class IntParamType>
    void fit(
        ValueType beta,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        GType& g,
        const WType& w,
        IntType ne,
        IntType nx,
        const XType& x,
        IntType nlam,
        ValueType flmin,
        const ULamType& ulam,
        ValueType thr,
        IntType maxit,
        const XMType& xm,
        const XSType& xs,
        const XVType& xv,
        IntType& lmu,
        AOType& ao,
        IAType& ia,
        KinType& kin,
        RSQOType& rsqo,
        ALMOType& almo,
        IntType& nlp,
        IntType& jerr,
        SetpbFType setpb_f,
        IntParamType int_param) const
    {
        FitPack<
            ValueType
            , JUType
            , VPType
            , CLType
            , GType
            , WType
            , IntType
            , XType
            , ULamType
            , XMType
            , XSType
            , XVType
            , AOType
            , IAType
            , KinType
            , RSQOType
            , ALMOType
            , SetpbFType
            , IntParamType> pack
        (beta, ju, vp, cl, ne, nx, x, nlam, flmin,
         ulam, thr, maxit, xv, lmu, ao, ia, kin, rsqo, almo, nlp, jerr, setpb_f, int_param, g,
         w, xm, xs);
        base_t::fit(pack);
    }

    template <class PackType>
    auto get_elnet_point(PackType&& pack) const 
    {
        return elnet_point_t(
                pack.thr, pack.maxit, pack.nx,
                pack.nlp, pack.ia, pack.g, pack.w, pack.x, 
                pack.xm, pack.xs, pack.xv, pack.vp, pack.cl, pack.ju);
    }

    template <class IntType
            , class PackType
            , class PathConfigPackType
            , class ElnetPointType>
    auto initialize_point(
            IntType m, 
            PackType&& pack, 
            PathConfigPackType&& path_pack,
            ElnetPointType&&) const
    {
        return gaussian_base_t::initialize_point(m, pack, path_pack, pack.g);    
    }
};

} // namespace glmnetpp
