#pragma once
#include <glmnetpp_bits/util/type_traits.hpp>
#include <glmnetpp_bits/elnet_driver/decl.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_path/binomial_multi_class.hpp>
#include <glmnetpp_bits/elnet_path/binomial_multi_class_group.hpp>
#include <glmnetpp_bits/elnet_path/sp_binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_path/sp_binomial_multi_class.hpp>
#include <glmnetpp_bits/elnet_path/sp_binomial_multi_class_group.hpp>
#include <glmnetpp_bits/util/types.hpp>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <glmnetpp_bits/elnet_driver/base.hpp>
#include <Eigen/Core>
#include <vector>

namespace glmnetpp {
namespace details {

template <bool is_dense>
struct FitPathBinomial
{
    template <class ValueType
            , class XType
            , class YType
            , class GType
            , class WType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class ULamType
            , class XMType
            , class XSType
            , class XVType
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
            const YType& y,
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
            const XVType& xv,
            ValueType thr,
            bool isd,
            bool intr,
            IntType maxit,
            IntType kopt,
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
        constexpr util::glm_type glm = util::glm_type::binomial;
        using mode_t = util::mode_type<glm>;

        if (g.cols() == 1) {
            auto y_1 = y.col(0);
            auto g_1 = g.col(0);
            Eigen::Map<Eigen::MatrixXd> ca_slice(ca.data(), nx, nlam);
            Eigen::Map<Eigen::VectorXd> a0_slice(a0.data(), a0.size());
            ElnetPath<glm, mode_t::two_class> elnet_path;
            elnet_path.fit(parm,ju,vp,cl,ne,nx,x,y_1,g_1,w,nlam,flmin,ulam,
                 thr,isd,intr,maxit,kopt,lmu,a0_slice,ca_slice,ia,nin,dev0,dev,alm,nlp,jerr,setpb_f, int_param);
        } 
        else if (kopt == 2) {
            ElnetPath<glm, mode_t::multi_class_group> elnet_path;
            elnet_path.fit(parm,ju,vp,cl,ne,nx,x,y,g,w,nlam,flmin,ulam,
                 thr,intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,setpb_f, int_param);
        }
        else {
            ElnetPath<glm, mode_t::multi_class> elnet_path;
            elnet_path.fit(parm,ju,vp,cl,ne,nx,x,y,g,w,nlam,flmin,ulam,thr,
                 isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,setpb_f,int_param);
        }
    }
};

template <>
struct FitPathBinomial<false>
{
    template <class ValueType
            , class XType
            , class YType
            , class GType
            , class WType
            , class JUType
            , class VPType
            , class CLType
            , class IntType
            , class ULamType
            , class XMType
            , class XSType
            , class XVType
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
            const YType& y,
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
            const XVType& xv,
            ValueType thr,
            bool isd,
            bool intr,
            IntType maxit,
            IntType kopt,
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
        constexpr util::glm_type glm = util::glm_type::binomial;
        using mode_t = util::mode_type<glm>;

        if (g.cols() == 1) {
            auto y_1 = y.col(0);
            auto g_1 = g.col(0);
            Eigen::Map<Eigen::MatrixXd> ca_slice(ca.data(), nx, nlam);
            Eigen::Map<Eigen::VectorXd> a0_slice(a0.data(), a0.size());
            SpElnetPath<glm, mode_t::two_class> elnet_path;
            elnet_path.fit(parm,ju,vp,cl,ne,nx,x,y_1,g_1,w,nlam,flmin,ulam,xm,xs,
                 thr,isd,intr,maxit,kopt,lmu,a0_slice,ca_slice,ia,nin,dev0,dev,alm,nlp,jerr,setpb_f, int_param);
        } 
        else if (kopt == 2) {
            SpElnetPath<glm, mode_t::multi_class_group> elnet_path;
            elnet_path.fit(parm,ju,vp,cl,ne,nx,x,y,g,w,nlam,flmin,ulam,
                 thr,intr,maxit,xm,xs,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,setpb_f,int_param);
        }
        else {
            SpElnetPath<glm, mode_t::multi_class> elnet_path;
            elnet_path.fit(parm,ju,vp,cl,ne,nx,x,y,g,w,nlam,flmin,ulam,xm,xs,thr,
                 isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,setpb_f,int_param);
        }
    }
};

} // namespace details

template <>
struct ElnetDriver<util::glm_type::binomial>
    : ElnetDriverBase
{
private:
    static constexpr util::glm_type glm = util::glm_type::binomial;
    using mode_t = util::mode_type<glm>;

public:
    template <class ValueType
            , class XType
            , class YType
            , class GType
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
            YType& y,
            GType& g,
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
            IntType kopt,
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
            IntParamType int_param
            ) const
    {
        using value_t = ValueType;
        using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
        using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;
        using bvec_t = std::vector<bool>;
        using x_t = std::decay_t<XType>;
        constexpr bool do_dense = util::is_dense<x_t>::value;
        using chkvars_t = typename std::conditional<do_dense,
                            Chkvars, SpChkvars>::type;
        using standardize_t = typename std::conditional<do_dense,
                            LStandardize1, SpLStandardize2>::type;
        using standardize_group_t = typename std::conditional<do_dense,
                            MultLStandardize1, MultSpLStandardize2>::type;

        try {
            vec_t vq = vp;
            this->normalize_penalty(vq);

            auto no = x.rows();
            auto ni = x.cols();
            auto nc = g.cols();

            vec_t ww(no);
            vec_t xm(ni); xm.setZero();
            vec_t xv; // only used if kopt == 2
            vec_t xs; // only used if isd > 0
            bvec_t ju(ni, false);

            if (kopt == 2) { xv.resize(ni); }

            if (do_dense) { 
                if (isd) { xs.setZero(ni); }
            } else {
                xs.setZero(ni);
            }

            chkvars_t::eval(x, ju);
            this->init_inclusion(jd, ju);

            for (int i = 0; i < no; ++i) {
                ww(i) = y.row(i).sum();
                if (ww(i)) y.row(i) /= ww(i);
            }
            auto sw = ww.sum();
            ww /= sw;

            if (nc != 1 && kopt == 2) {
                standardize_group_t::eval(x, ww, ju, isd, intr, xm, xs, xv);
            }
            else {
                standardize_t::eval(x, ww, ju, isd, intr, xm, xs);
            }
            if (isd) { 
                for (int j = 0; j < ni; ++j) cl.col(j) *= xs(j);
            }

            details::FitPathBinomial<do_dense>::eval(
                    parm, x, y, g, ww, ju, vq, cl, ne, nx,
                    nlam, flmin, ulam, xm, xs, xv, thr, isd, intr, maxit, kopt,
                    lmu, a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, setpb_f, int_param
                    );

            if (jerr > 0) return; 

            dev0 *= 2.0*sw;
            for (int k = 0; k < lmu; ++k) {
                auto nk = nin(k);
                for (int ic = 0; ic < nc; ++ic) {
                    Eigen::Map<mat_t> ca_slice(ca.data() + k*nx*nc, nx, nc);
                    if (isd) { 
                        for (int l = 0; l < nk; ++l) {
                            ca_slice(l, ic) /= xs(ia(l)-1);
                        }
                    }
                    if (intr == 0) { a0(ic,k)=0.0; }
                    else { 
                        for (int i = 0; i < nk; ++i) {
                            a0(ic, k) -= ca_slice(i, ic) * xm(ia(i)-1);
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
