#pragma once
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/binomial_base.hpp>

namespace glmnetpp {

template <class ElnetPointPolicy>
struct ElnetPath<
    util::glm_type::binomial,
    util::mode_type<util::glm_type::binomial>::multi_class,
    ElnetPointPolicy>
        : ElnetPathBinomialMultiClassBase
        , ElnetPathCRTPBase<
            ElnetPath<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::multi_class,
                ElnetPointPolicy> >
{
private:
    using base_t = ElnetPathBinomialMultiClassBase;
    using crtp_base_t = ElnetPathCRTPBase<
            ElnetPath<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::multi_class,
                ElnetPointPolicy> >;
    using elnet_point_t = ElnetPointPolicy;

public:
    using base_t::initialize_path;
    using base_t::process_point_fit;
    using base_t::process_path_fit;

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
    void fit(
        ValueType beta,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        IntType ne,
        IntType nx,
        const XType& x,
        const YType& y,
        GType& g,
        const WType& w,
        IntType nlam,
        ValueType flmin,
        const ULamType& ulam,
        ValueType thr,
        bool isd,
        bool intr,
        IntType maxit,
        IntType kopt,
        IntType& lmu,
        A0Type& a0,
        AOType& ao,
        IAType& ia,
        KinType& kin,
        ValueType& dev0,
        DevType& dev,
        ALMType& alm,
        IntType& nlp,
        IntType& jerr,
        SetpbFType setpb_f,
        const IntParamType& int_param) const
    {
        FitPack<
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
            , IntParamType> pack{
                // build sub-pack
                {beta, ju, vp, cl, ne, nx, x, nlam, flmin,
                 ulam, thr, maxit, lmu, ao, ia, kin, alm, nlp, jerr, setpb_f, int_param},
                // add new members
                y, g, w, isd, intr, kopt, a0,  dev0, dev
        };
        crtp_base_t::fit(pack);
    }

    template <class FitPackType, class PathConfigPackType>
    elnet_point_t get_elnet_point(
            const FitPackType& pack, 
            PathConfigPackType&& path_pack) const 
    {
        auto& sp = pack.sub_pack;
        return elnet_point_t(
                pack.isd, pack.intr, pack.kopt, sp.thr, sp.maxit,
                sp.nx, sp.nlp, sp.ia, pack.g, pack.dev0, 
                sp.x, pack.y, pack.w, sp.vp, sp.cl, sp.ju, 
                path_pack.is, sp.int_param);
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
        return base_t::initialize_point(m, lmda_curr, pack, path_pack, elnet_point.abs_grad());
    }

};

} // namespace glmnetpp
