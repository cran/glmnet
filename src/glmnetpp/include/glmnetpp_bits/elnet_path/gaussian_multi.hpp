#pragma once 
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_base.hpp>

namespace glmnetpp {

/* 
 * Gaussian multi-response method elastic net path-solver.
 */
template <class ElnetPointPolicy>
struct ElnetPath<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::multi,
    ElnetPointPolicy>
        : ElnetPathGaussianMultiBase
        , ElnetPathCRTPBase<
            ElnetPath<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::multi,
                ElnetPointPolicy> >
{
private:
    using base_t = ElnetPathGaussianMultiBase;
    using crtp_base_t = ElnetPathCRTPBase<
            ElnetPath<util::glm_type::gaussian, 
                      util::mode_type<util::glm_type::gaussian>::multi,
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
        ValueType ys0,
        IntType& lmu,
        AOType& ao,
        IAType& ia,
        KinType& kin,
        RSQOType& rsqo,
        ALMOType& almo,
        IntType& nlp,
        IntType& jerr,
        SetpbFType setpb_f,
        const IntParamType& int_param) const
    {
        base_t::FitPack<
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
        {
            // build sub-pack
            {
                // build sub-pack
                {beta, ju, vp, cl, ne, nx, x, nlam, flmin,
                 ulam, thr, maxit, lmu, ao, ia, kin, almo, nlp, jerr, setpb_f, int_param},
                // add new members
                xv, rsqo
            }, 
            // add new members
            y, ys0
        };
        crtp_base_t::fit(pack);
    }

    template <class FitPackType, class PathConfigPackType>
    auto get_elnet_point(const FitPackType& pack, const PathConfigPackType&) const 
    {
        auto& sp = pack.sub_pack;
        auto& ssp = sp.sub_pack;
        return elnet_point_t(
                ssp.thr, ssp.maxit, ssp.nx, ssp.nlp, ssp.ia, pack.ys0, pack.y, ssp.x, 
                sp.xv, ssp.vp, ssp.cl, ssp.ju, ssp.int_param);
    }
};

} // namespace glmnetpp


