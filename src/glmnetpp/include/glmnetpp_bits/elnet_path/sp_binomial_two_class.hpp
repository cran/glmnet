#pragma once 
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/binomial_base.hpp>
#include <glmnetpp_bits/util/macros.hpp>

namespace glmnetpp {

template <class ElnetPointPolicy>
struct SpElnetPath<
    util::glm_type::binomial,
    util::mode_type<util::glm_type::binomial>::two_class,
    ElnetPointPolicy>
        : ElnetPathBinomialTwoClassBase
        , ElnetPathCRTPBase<
            SpElnetPath<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::two_class,
                ElnetPointPolicy> >
{
private:
    using base_t = ElnetPathBinomialTwoClassBase;
    using crtp_base_t = ElnetPathCRTPBase<
                SpElnetPath<
                      util::glm_type::binomial, 
                      util::mode_type<util::glm_type::binomial>::two_class,
                      ElnetPointPolicy> >;
    using elnet_point_t = ElnetPointPolicy;
    using typename base_t::state_t;

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
            , class XBType
            , class XSType
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
        const XBType& xb;
        const XSType& xs;
    };

public:

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
            , class XBType
            , class XSType
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
        const XBType& xb,
        const XSType& xs,
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
            , XBType
            , XSType
            , A0Type
            , AOType
            , IAType
            , KinType
            , DevType
            , ALMType
            , SetpbFType
            , IntParamType> pack{
            // build sub-pack
            {
                // build sub-pack
                {
                    beta, ju, vp, cl, ne, nx, x, nlam, flmin,
                    ulam, thr, maxit, lmu, ao, ia, kin, alm, nlp, jerr, setpb_f, int_param
                },
                // add new members
                y, g, w, isd, intr, kopt, a0,  dev0, dev
            },
            xb, xs
        };
        crtp_base_t::fit(pack);
    }

    template <class FitPackType, class PathConfigPackType>
    elnet_point_t get_elnet_point(
            const FitPackType& pack, 
            const PathConfigPackType&) const 
    {
        auto& sp = pack.sub_pack;
        auto& ssp = sp.sub_pack;
        return elnet_point_t(
                sp.isd, sp.intr, sp.kopt, ssp.thr, ssp.maxit,
                ssp.nx, ssp.nlp, ssp.ia, sp.g, sp.dev0, 
                ssp.x, sp.y, sp.w, pack.xb, pack.xs, ssp.vp, ssp.cl, ssp.ju, 
                ssp.int_param);
    }

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
    auto initialize_point(
            IntType m, 
            ValueType&& lmda_curr,
            const FitPackType& pack, 
            const PathConfigPackType& path_pack,
            const ElnetPointType& elnet_point) const
    {
        return base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, elnet_point);
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
        return base_t::process_point_fit(pack.sub_pack, path_pack, point_pack, elnet_point);
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
