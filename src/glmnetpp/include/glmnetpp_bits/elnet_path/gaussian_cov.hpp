#pragma once 
#include <glmnetpp_bits/elnet_path/decl.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_base.hpp>

namespace glmnetpp {

/* 
 * Gaussian covariance method elastic net path-solver.
 */
template <class ElnetPointPolicy>
struct ElnetPath<
    util::glm_type::gaussian,
    util::mode_type<util::glm_type::gaussian>::cov,
    ElnetPointPolicy>
        : ElnetPathGaussianBase
        , ElnetPathCRTPBase<
            ElnetPath<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::cov,
                ElnetPointPolicy> >
{
private:
    using base_t = ElnetPathGaussianBase;
    using crtp_base_t = ElnetPathCRTPBase<
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
    struct FitPack
    {
        using sub_pack_t = typename base_t::template FitPack<
            ValueType, JUType, VPType, CLType, IntType,
            XType, ULamType, XVType, AOType, IAType, KinType,
            RSQOType, ALMOType, SetpbFType, IntParamType>;
        using value_t = typename sub_pack_t::value_t;
        using int_t = typename sub_pack_t::int_t;

        int_t& err_code() const { return sub_pack.err_code(); }
        int_t path_size() const { return sub_pack.path_size(); }

        sub_pack_t sub_pack;
        GType& g;
    };

public:
    using base_t::process_path_fit;

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
        const IntParamType& int_param) const
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
            g
        };
        crtp_base_t::fit(pack);
    }

    /*
     * Builds a point-solver using the arguments.
     *
     * @param   pack        object of FitPack of the current class.
     */
    template <class FitPackType, class PathConfigPackType>
    auto get_elnet_point(const FitPackType& pack, const PathConfigPackType&) const 
    {
        auto& sp = pack.sub_pack;
        auto& ssp = sp.sub_pack;
        return elnet_point_t(
                ssp.thr, ssp.maxit, ssp.nx, ssp.nlp, ssp.ia, pack.g, ssp.x, 
                sp.xv, ssp.vp, ssp.cl, ssp.ju);
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
            const ElnetPointType&) const
    {
        return base_t::initialize_point(m, lmda_curr, pack.sub_pack, path_pack, pack.g);
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
};

} // namespace glmnetpp
