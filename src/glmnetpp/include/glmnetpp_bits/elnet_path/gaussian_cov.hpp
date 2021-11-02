#pragma once 
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_base.hpp>

namespace glmnetpp {

template <class ElnetPointPolicy>
struct ElnetPath<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::cov,
    ElnetPointPolicy>
        : private ElnetPathGaussianBase
        , ElnetPathBase<
            ElnetPath<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::cov,
                ElnetPointPolicy> >
{
private:
    using gaussian_base_t = ElnetPathGaussianBase;
    using base_t = ElnetPathBase<
            ElnetPath<util::glm_type::gaussian, 
                      util::mode_type<util::glm_type::gaussian>::cov,
                      ElnetPointPolicy> >;
    using elnet_point_t = ElnetPointPolicy;

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class GType
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
            GType& _g
               )
            : base_t{beta, ju, vp, cl, ne, nx, x, nlam, flmin, ulam, thr, maxit, xv, 
                     lmu, ao, ia, kin, rsqo, almo, nlp, jerr, setpb_f, int_param}
            , g(_g)
        {}
    };

public:
    using gaussian_base_t::initialize_path;

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class GType
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
    void fit(
        ValueType beta,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        GType& g,
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
        IntParamType int_param) const
    {
        FitPack<
            ValueType
            , JUType
            , VPType
            , CLType
            , GType
            , IntType
            , XType
            , ULamType
            , XVType
            , AOType
            , IAType
            , KinType
            , RSQOType
            , ALMOType
            , SetpbFType
            , IntParamType> pack
        (beta, ju, vp, cl, ne, nx, x, nlam, flmin,
         ulam, thr, maxit, xv, lmu, ao, ia, kin, rsqo, almo, nlp, jerr, setpb_f, int_param, g);
        base_t::fit(pack);
    }

    template <class PackType>
    auto get_elnet_point(PackType&& pack) const 
    {
        return elnet_point_t(
                pack.thr, pack.maxit, pack.nx,
                pack.nlp, pack.ia, pack.g, pack.x, 
                pack.xv, pack.vp, pack.cl, pack.ju);
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
