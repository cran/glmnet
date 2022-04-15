#pragma once 
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/poisson_base.hpp>
#include <glmnetpp_bits/util/macros.hpp>

namespace glmnetpp {

template <class ElnetPointPolicy>
struct ElnetPath<
    util::glm_type::poisson,
    util::mode_type<util::glm_type::poisson>::naive,
    ElnetPointPolicy>
        : ElnetPathPoissonBase
        , ElnetPathCRTPBase<
            ElnetPath<
                util::glm_type::poisson,
                util::mode_type<util::glm_type::poisson>::naive,
                ElnetPointPolicy> >
{
private:
    using base_t = ElnetPathPoissonBase;
    using crtp_base_t = ElnetPathCRTPBase<
            ElnetPath<util::glm_type::poisson, 
                      util::mode_type<util::glm_type::poisson>::naive,
                      ElnetPointPolicy> >;
    using elnet_point_t = ElnetPointPolicy;

public:
    using base_t::initialize_path;
    using base_t::initialize_point;
    using base_t::process_path_fit;
    using base_t::process_point_fit;

    template <class ValueType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class XType
            , class YType
            , class GType
            , class QType
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
        const QType& q,
        IntType nlam,
        ValueType flmin,
        const ULamType& ulam,
        ValueType thr,
        bool intr,
        IntType maxit,
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
            , QType
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
                y, g, q, intr, a0, dev0, dev
        };
        crtp_base_t::fit(pack);
    }

    template <class FitPackType, class PathConfigPackType>
    elnet_point_t get_elnet_point(
            const FitPackType& pack, 
            const PathConfigPackType&) const 
    {
        auto& sp = pack.sub_pack;
        return elnet_point_t(
                pack.intr, sp.thr, sp.maxit,
                sp.nx, sp.nlp, sp.ia, pack.dev0, 
                sp.x, pack.y, pack.g, pack.q, sp.vp, sp.cl, sp.ju, 
                sp.int_param);
    }
};

} // namespace glmnetpp
