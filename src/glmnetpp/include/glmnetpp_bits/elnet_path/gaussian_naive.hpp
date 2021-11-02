#pragma once 
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_base.hpp>

namespace glmnetpp {

template <class ElnetPointPolicy>
struct ElnetPath<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::naive,
    ElnetPointPolicy>
        : private ElnetPathGaussianBase
        , ElnetPathBase<
            ElnetPath<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::naive,
                ElnetPointPolicy> >
{
private:
    using gaussian_base_t = ElnetPathGaussianBase;
    using base_t = ElnetPathBase<
            ElnetPath<util::glm_type::gaussian, 
                      util::mode_type<util::glm_type::gaussian>::naive,
                      ElnetPointPolicy> >;
    using elnet_point_t = ElnetPointPolicy;

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
        YType& y;

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
            YType& _y
               )
            : base_t{beta, ju, vp, cl, ne, nx, x, nlam, flmin, ulam, thr, maxit, xv, 
                     lmu, ao, ia, kin, rsqo, almo, nlp, jerr, setpb_f, int_param}
            , y(_y)
        {}
    };

public:
    using gaussian_base_t::initialize_path;

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
    void fit(
        ValueType beta,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        YType& y,
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
            , YType
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
         ulam, thr, maxit, xv, lmu, ao, ia, kin, rsqo, almo, nlp, jerr, setpb_f, int_param, y);
        base_t::fit(pack);
    }

    template <class PackType>
    elnet_point_t get_elnet_point(PackType&& pack) const 
    {
        return elnet_point_t(
                pack.thr, 
                pack.maxit, 
                pack.nx,
                pack.nlp, 
                pack.ia, 
                pack.y, 
                pack.x, 
                pack.xv, 
                pack.vp, 
                pack.cl, 
                pack.ju);
    }

    template <class IntType
            , class PackType
            , class PathConfigPackType
            , class ElnetPointType>
    auto initialize_point(
            IntType m, 
            PackType&& pack, 
            PathConfigPackType&& path_pack,
            ElnetPointType&& elnet_point) const
    {
        return gaussian_base_t::initialize_point(m, pack, path_pack, elnet_point.abs_grad());
    }

};

} // namespace glmnetpp


