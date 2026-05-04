#pragma once
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/elnet_driver/decl.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/cox_naive.hpp>
#include <glmnetpp_bits/elnet_path/sp_cox_naive.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <glmnetpp_bits/elnet_driver/base.hpp>
#include <Eigen/Core>
#include <vector>

namespace glmnetpp {
namespace details {

template <bool is_dense>
struct FitPathCox
{
    template <class ValueType
            , class XType
            , class SurvType      // Survival response (event times, status, etc.)
            , class GType
            , class WType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class ULamType
            , class XMType
            , class XSType
            , class LmuType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class DevType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    static void eval(
            ValueType parm,
            const XType& x,
            const SurvType& surv,   // Survival data instead of y
            GType& g,               // Linear predictor (eta)
            const WType& w,
            const JUType& ju,
            const VPType& vp,
            const CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const ULamType& ulam,
            const XMType&,
            const XSType&,
            ValueType thr,
            IntType maxit,
            LmuType& lmu,
            A0Type& a0,             // Not used for Cox (no intercept) but kept for interface
            CAType& ca,
            IAType& ia,
            NinType& nin,
            ValueType& dev0,
            DevType& dev,
            AlmType& alm,
            IntType& nlp,
            IntType& jerr,
            SetpbFType setpb_f,
            const IntParamType& int_param
            )
    {
        constexpr util::glm_type glm = util::glm_type::cox;
        using mode_t = util::mode_type<glm>;

        ElnetPath<glm, mode_t::naive> elnet_path;
        elnet_path.fit(parm, ju, vp, cl, ne, nx, x, surv, g, w, nlam, flmin, ulam,
             thr, maxit, lmu, a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, setpb_f, int_param);
    }
};

// Sparse specialization
template <>
struct FitPathCox<false>
{
    template <class ValueType
            , class XType
            , class SurvType
            , class GType
            , class WType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class ULamType
            , class XMType
            , class XSType
            , class LmuType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class DevType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    static void eval(
            ValueType parm,
            const XType& x,
            const SurvType& surv,
            GType& g,
            const WType& w,
            const JUType& ju,
            const VPType& vp,
            const CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const ULamType& ulam,
            const XMType& xm,
            const XSType& xs,
            ValueType thr,
            IntType maxit,
            LmuType& lmu,
            A0Type& a0,
            CAType& ca,
            IAType& ia,
            NinType& nin,
            ValueType& dev0,
            DevType& dev,
            AlmType& alm,
            IntType& nlp,
            IntType& jerr,
            SetpbFType setpb_f,
            const IntParamType& int_param
            )
    {
        constexpr util::glm_type glm = util::glm_type::cox;
        using mode_t = util::mode_type<glm>;

        SpElnetPath<glm, mode_t::naive> elnet_path;
        elnet_path.fit(parm, ju, vp, cl, ne, nx, x, surv, g, w, nlam, flmin, ulam,
             xm, xs, thr, maxit, lmu, a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, setpb_f, int_param);
    }
};

} // namespace details

template <>
struct ElnetDriver<util::glm_type::cox>
    : ElnetDriverBase
{
private:
    static constexpr util::glm_type glm = util::glm_type::cox;
    using mode_t = util::mode_type<glm>;

public:
    template <class ValueType
            , class XType
            , class SurvType      // Survival response type
            , class GType
            , class WType
            , class JDType
            , class VPType
            , class CLType
            , class ULamType
            , class IntType
            , class LmuType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class DevType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    void fit(
            ValueType parm,
            XType& x,
            SurvType& surv,       // Survival data
            GType& g,             // Linear predictor (eta)
            const WType& w,
            const JDType& jd,
            const VPType& vp,
            CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const ULamType& ulam,
            ValueType thr,
            bool isd,             // Standardize x?
            IntType maxit,
            LmuType& lmu,
            A0Type& a0,           // Intercepts (not used for Cox, but kept for consistency)
            CAType& ca,
            IAType& ia,
            NinType& nin,
            ValueType& dev0,
            DevType& dev,
            AlmType& alm,
            IntType& nlp,
            IntType& jerr,
            SetpbFType setpb_f,
            IntParamType int_param
            ) const
    {
        using value_t = ValueType;
        using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
        using bvec_t = std::vector<bool>;
        using x_t = std::decay_t<XType>;
        constexpr bool do_dense = util::is_dense<x_t>::value;
        using chkvars_t = typename std::conditional<do_dense,
                            Chkvars, SpChkvars>::type;
        using standardize_t = typename std::conditional<do_dense,
                            LStandardize1, SpLStandardize2>::type;

        try {
            vec_t vq = vp;
            this->normalize_penalty(vq);

            auto no = x.rows();
            auto ni = x.cols();

            vec_t ww(no);
            vec_t xm(ni); xm.setZero();
            vec_t xs;   // only used if isd is true
            bvec_t ju(ni, false);

            if (do_dense) {
                if (isd) { xs.setZero(ni); }
            } else {
                xs.setZero(ni);
            }

            chkvars_t::eval(x, ju);
            this->init_inclusion(jd, ju);

            ww.array() = w.array().max(0.0);
            auto sw = ww.sum();
            if (sw <= 0.0) throw util::non_positive_weight_sum_error();
            ww /= sw;

            // Standardize x (Cox has no intercept, so intr=false)
            constexpr bool intr = false;
            standardize_t::eval(x, ww, ju, isd, intr, xm, xs);
            if (isd) {
                for (int j = 0; j < ni; ++j) {
                    cl.col(j) *= xs(j);
                }
            }

            details::FitPathCox<do_dense>::eval(
                    parm, x, surv, g, ww, ju, vq, cl, ne, nx,
                    nlam, flmin, ulam, xm, xs, thr, maxit,
                    lmu, a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, setpb_f, int_param
                    );

            if (jerr > 0) return;

            // Note: Unlike other families, Cox deviance from coxdev is already in standard form
            // (-2 * log partial likelihood), so we don't multiply by 2.
            // We only scale by sum of weights for consistency with R output.
            dev0 *= sw;
            for (int k = 0; k < lmu; ++k) {
                auto nk = nin(k);
                if (isd) {
                    for (int l = 0; l < nk; ++l) {
                        ca(l, k) /= xs(ia(l)-1);
                    }
                }
                // Cox has no intercept, so a0 stays at 0
                a0(k) = 0.0;
            }
        }
        catch (const util::elnet_error& e) {
            jerr = e.err_code(0);
        }
    }
};

} // namespace glmnetpp
