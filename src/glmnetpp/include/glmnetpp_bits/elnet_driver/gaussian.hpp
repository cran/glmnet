#pragma once
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/elnet_driver/decl.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/chkvars.hpp>
#include <glmnetpp_bits/standardize.hpp>
#include <Eigen/Core>
#include <vector>

namespace glmnetpp {
namespace details {

template <bool is_dense>
struct FitPath
{
    template <class ValueType
            , class XType
            , class YType
            , class GType
            , class WType
            , class JUType
            , class VQType
            , class XMType
            , class XSType
            , class XVType
            , class CLType
            , class VLamType
            , class IntType
            , class LmuType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class RsqType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    static void eval(
            bool ka,
            ValueType parm,
            const XType& x,
            YType& y,
            GType& g,
            const WType& w,
            const JUType& ju,
            const VQType& vq,
            const XMType&,
            const XSType&,
            const XVType& xv,
            const CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const VLamType& vlam,
            ValueType thr,
            bool isd,
            bool intr,
            IntType maxit,
            LmuType& lmu,
            A0Type& a0,
            CAType& ca,
            IAType& ia,
            NinType& nin,
            RsqType& rsq,
            AlmType& alm,
            IntType& nlp,
            IntType& jerr,
            SetpbFType setpb_f,
            IntParamType int_param
            ) 
    {
        constexpr util::glm_type glm = util::glm_type::gaussian;
        using mode_t = util::mode_type<glm>;

        // cov method
        if (!ka) {
            ElnetPath<glm, mode_t::cov> elnet_path;
            elnet_path.fit(
                    parm, ju, vq, cl, g, ne, nx, x, nlam, flmin, vlam, thr, maxit, xv,
                    lmu, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
        } 
        // naive method
        else {
            ElnetPath<glm, mode_t::naive> elnet_path;
            elnet_path.fit(
                    parm, ju, vq, cl, y, ne, nx, x, nlam, flmin, vlam, thr, maxit, xv,
                    lmu, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
        }
    }
};

template <>
struct FitPath<false>
{
    template <class ValueType
            , class XType
            , class YType
            , class GType
            , class WType
            , class JUType
            , class VQType
            , class XMType
            , class XSType
            , class XVType
            , class CLType
            , class VLamType
            , class IntType
            , class LmuType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class RsqType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    static void eval(
            bool ka,
            ValueType parm,
            const XType& x,
            YType& y,
            GType& g,
            const WType& w,
            const JUType& ju,
            const VQType& vq,
            const XMType& xm,
            const XSType& xs,
            const XVType& xv,
            const CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const VLamType& vlam,
            ValueType thr,
            bool isd,
            bool intr,
            IntType maxit,
            LmuType& lmu,
            A0Type& a0,
            CAType& ca,
            IAType& ia,
            NinType& nin,
            RsqType& rsq,
            AlmType& alm,
            IntType& nlp,
            IntType& jerr,
            SetpbFType setpb_f,
            IntParamType int_param
            ) 
    {
        constexpr util::glm_type glm = util::glm_type::gaussian;
        using mode_t = util::mode_type<glm>;

        // cov method
        if (!ka) {
            SpElnetPath<glm, mode_t::cov> elnet_path;
            elnet_path.fit(
                    parm, ju, vq, cl, g, w, ne, nx, x, nlam, flmin, vlam, thr, maxit, xm, xs, xv,
                    lmu, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
        } 
        // naive method
        else {
            SpElnetPath<glm, mode_t::naive> elnet_path;
            elnet_path.fit(
                    parm, ju, vq, cl, y, w, ne, nx, x, nlam, flmin, vlam, thr, maxit, xm, xs, xv,
                    lmu, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
        }
    }
};

} // namespace details

template <>
struct ElnetDriver<util::glm_type::gaussian>
{
private:
    static constexpr util::glm_type glm = util::glm_type::gaussian;
    using mode_t = util::mode_type<glm>;

public:
    template <class ValueType
            , class XType
            , class YType
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
            , class RsqType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    void fit(
            bool ka,
            ValueType parm,
            XType& x,
            YType& y,
            WType& w,
            const JDType& jd,
            const VPType& vp,
            CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const ULamType& ulam,
            ValueType thr,
            bool isd,
            bool intr,
            IntType maxit,
            LmuType& lmu,
            A0Type& a0,
            CAType& ca,
            IAType& ia,
            NinType& nin,
            RsqType& rsq,
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
        using standardize_cov_t = typename std::conditional<do_dense,
                            Standardize, SpStandardize>::type;
        using standardize_naive_t = typename std::conditional<do_dense,
                            Standardize1, SpStandardize1>::type;

        if (vp.maxCoeff() <= 0) { jerr = 10000; return; }

        auto ni = x.cols();
        vec_t vq = vp.unaryExpr([](auto x){ return std::max(x, 0.); }); 
        vq *= ni / vq.sum();

        vec_t g;            // only naive version uses it
        vec_t xm(ni); xm.setZero();
        vec_t xs(ni); xs.setZero();
        vec_t xv(ni); xv.setZero();
        vec_t vlam(nlam); vlam.setZero();
        bvec_t ju(ni, false);

        chkvars_t::eval(x, ju);

        if (jd(0) > 0) {
            for (int i = 1; i < jd(0) + 1; ++i) {
                ju[jd(i)-1] = false;
            }
        }

        // can't find true value in ju
        if (std::find_if(ju.begin(), ju.end(), [](auto x) { return x;}) == ju.end()) {
            jerr=7777; 
            return;
        } 

        value_t ym = 0; 
        value_t ys = 0;

        // cov method
        if (!ka) {
            g.setZero(ni);
            standardize_cov_t::eval(x, y, w, isd, intr, ju, g, xm, xs, ym, ys, xv);
        } 

        // naive method
        else {
            standardize_naive_t::eval(x, y, w, isd, intr, ju, xm, xs, ym, ys, xv);
        }

        cl /= ys; 
        if (isd) { 
            for (int j = 0; j < ni; ++j) {
                cl.col(j) *= xs(j);
            }
        }

        if (flmin >= 1.0) vlam = ulam / ys;

        details::FitPath<do_dense>::eval(
                ka, parm, x, y, g, w, ju, vq, xm, xs, xv, cl, ne, nx,
                nlam, flmin, vlam, thr, isd, intr, maxit,
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param
                );
            
        if (jerr > 0) return;

        for (int k = 0; k < lmu; ++k) {
            alm(k) *= ys; 
            auto nk = nin(k);
            for (int l = 0; l < nk; ++l) {
                ca(l,k) *= ys / xs(ia(l)-1);
            }
            a0(k)=0.0;
            if (intr) {
                for (int i = 0; i < nk; ++i) {
                    a0(k) -= ca(i, k) * xm(ia(i)-1);
                }
                a0(k) += ym;
            }
        }
    }
};

} // namespace glmnetpp
