#pragma once
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/elnet_driver/decl.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_multi.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <glmnetpp_bits/elnet_driver/base.hpp>
#include <Eigen/Core>
#include <vector>

namespace glmnetpp {
namespace details {

/*
 * Primary definition: dense, not multi
 */
template <bool is_dense, bool is_multi>
struct FitPathGaussian
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

/*
 * Specialization: sparse, not multi
 */
template <>
struct FitPathGaussian<false, false>
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


/*
 * Specialization: dense, multi
 */
template <>
struct FitPathGaussian<true, true>
{
    template <class ValueType
            , class JUType
            , class VQType
            , class CLType
            , class YType
            , class WType
            , class IntType
            , class XType
            , class VLamType
            , class XMType
            , class XSType
            , class XVType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class RsqType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    static void eval(
            ValueType parm,
            const JUType& ju,
            const VQType& vq,
            const CLType& cl,
            YType& y,
            const WType&,
            IntType ne,
            IntType nx,
            const XType& x,
            IntType nlam,
            ValueType flmin,
            const VLamType& vlam,
            ValueType thr,
            IntType maxit,
            const XMType&,
            const XSType&,
            const XVType& xv,
            ValueType ys0,
            IntType& lmu,
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

        ElnetPath<glm, mode_t::multi> elnet_path;
        elnet_path.fit(
                parm, ju, vq, cl, y, ne, nx, x, nlam, flmin, vlam, thr, maxit, xv,
                ys0, lmu, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
    }
};

/*
 * Specialization: sparse, multi
 */
template <>
struct FitPathGaussian<false, true>
{
    template <class ValueType
            , class JUType
            , class VQType
            , class CLType
            , class YType
            , class WType
            , class IntType
            , class XType
            , class VLamType
            , class XMType
            , class XSType
            , class XVType
            , class A0Type
            , class CAType
            , class IAType
            , class NinType
            , class RsqType
            , class AlmType
            , class SetpbFType
            , class IntParamType>
    static void eval(
            ValueType parm,
            const JUType& ju,
            const VQType& vq,
            const CLType& cl,
            YType& y,
            const WType& w,
            IntType ne,
            IntType nx,
            const XType& x,
            IntType nlam,
            ValueType flmin,
            const VLamType& vlam,
            ValueType thr,
            IntType maxit,
            const XMType& xm,
            const XSType& xs,
            const XVType& xv,
            ValueType ys0,
            IntType& lmu,
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

        SpElnetPath<glm, mode_t::multi> elnet_path;
        elnet_path.fit(
                parm, ju, vq, cl, y, w, ne, nx, x, nlam, flmin, vlam, thr, maxit, xm, xs, xv,
                ys0, lmu, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
    }
};

} // namespace details

template <>
struct ElnetDriver<util::glm_type::gaussian>
    : ElnetDriverBase
{
private:
    static constexpr util::glm_type glm = util::glm_type::gaussian;
    using mode_t = util::mode_type<glm>;

public:

    /* 
     * Fit function for dense/sparse naive/cov.
     */
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

        try {
            vec_t vq = vp;
            this->normalize_penalty(vq);

            auto ni = x.cols();

            vec_t g;            // only naive version uses it
            vec_t xm(ni); xm.setZero();
            vec_t xs(ni); xs.setZero();
            vec_t xv(ni); xv.setZero();
            vec_t vlam(nlam); vlam.setZero();
            bvec_t ju(ni, false);

            chkvars_t::eval(x, ju);
            this->init_inclusion(jd, ju);

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

            details::FitPathGaussian<do_dense, false>::eval(
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
        catch (const util::elnet_error& e) {
            jerr = e.err_code(0);
        }
    }

    /*
     * Fit function for dense/sparse multi-response.
     */
    template <class ValueType
            , class XType
            , class YType
            , class WType
            , class JDType
            , class VPType
            , class CLType
            , class IntType
            , class ULamType
            , class A0Type
            , class AOType
            , class IAType
            , class KinType
            , class RsqoType
            , class AlmoType
            , class SetpbFType
            , class IntParamType
            >
    void fit(
            ValueType beta,
            XType& x,
            YType& y,
            WType& w,
            const JDType& jd,
            const VPType& vp,
            const CLType& cl,
            IntType ne,
            IntType nx,
            IntType nlam,
            ValueType flmin,
            const ULamType& ulam,
            ValueType thr,
            bool isd, 
            bool jsd,
            bool intr,
            IntType maxit,
            IntType& lmu,
            A0Type& a0,
            AOType& ca,
            IAType& ia,
            KinType& nin,
            RsqoType& rsq,
            AlmoType& alm,
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
                            MultStandardize1, MultSpStandardize1>::type;

        try {
            vec_t vq = vp;
            this->normalize_penalty(vq);

            auto ni = x.cols();
            auto nr = y.cols();

            vec_t clt(2 * nr * ni);
            vec_t xm(ni);
            vec_t xs(ni);
            vec_t ym(nr);
            vec_t ys(nr);
            bvec_t ju(ni, false);
            vec_t xv(ni);

            chkvars_t::eval(x, ju);
            this->init_inclusion(jd, ju);

            value_t ys0 = 0.0;
            standardize_t::eval(x, y, w, isd, jsd, intr, ju, xm, xs, ym, ys, xv, ys0);

            for (int j = 0; j < ni; ++j) {
                Eigen::Map<Eigen::MatrixXd> clt_slice(
                        clt.data() + j * 2 * nr, 2, nr);
                for (int k = 0; k < nr; ++k) {
                    for (int i = 0; i < 2; ++i) {
                        clt_slice(i, k) = cl(i,j);
                        if (isd) clt_slice(i, k) *= xs(j);
                        if (jsd) clt_slice(i, k) /= ys(k);
                    }
                }
            }

            details::FitPathGaussian<do_dense, true>::eval(
                    beta, ju, vq, clt, y, w, ne, nx, x,
                    nlam, flmin, ulam, thr, maxit, xm, xs, xv, 
                    ys0, lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr, setpb_f, int_param);
                
            if (jerr > 0) return;

            for (int k = 0; k < lmu; ++k) {
                auto nk = nin(k);
                Eigen::Map<Eigen::MatrixXd> ca_slice(
                        ca.data() + k * nx * nr, nx, nr);
                for (int j = 0; j < nr; ++j) {
                    for (int l = 0; l < nk; ++l) {
                        ca_slice(l, j) *= ys(j) / xs(ia(l)-1);
                    }
                    if (!intr) a0(j,k) = 0.0;
                    else {
                        a0(j,k) = ym(j);
                        for (int l = 0; l < nk; ++l) {
                            a0(j,k) -= ca_slice(l,j) * xm(ia(l)-1);
                        }
                    }
                }
            }
        }
        catch (const util::elnet_error& e) {
            jerr = e.err_code(0);
        }
    }
};

} // namespace glmnetpp
