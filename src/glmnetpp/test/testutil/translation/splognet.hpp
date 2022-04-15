#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

// have been separately unittested
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>        
#include <glmnetpp_bits/elnet_driver/standardize.hpp> 
#include <testutil/translation/splognetn.hpp>
#include <testutil/translation/splognet2n.hpp>
#include <testutil/translation/multsplognetn.hpp>

namespace glmnetpp {
namespace transl {

template <class FloatType
        , bool apply_suggested_fixes
        , class ValueType
        , class XType
        , class YType
        , class GType
        , class JDType
        , class VPType
        , class CLType
        , class IntType
        , class ULamType
        , class A0Type
        , class CAType
        , class IAType
        , class NinType
        , class DevType
        , class AlmType
        >
inline void splognet(
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
        IntType& lmu,
        A0Type& a0,
        CAType& ca,
        IAType& ia,
        NinType& nin,
        ValueType& dev0,
        DevType& dev,
        AlmType& alm,
        IntType& nlp,
        IntType& jerr)
{
    using value_t = ValueType;
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using mat_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;
    using bvec_t = std::vector<bool>;

    // TODO: this error should be generalized. Gaussian too.
    if (vp.maxCoeff() <= 0) { jerr = 10000; return; }

    auto no = x.rows();
    auto ni = x.cols();
    auto nc = g.cols();
    vec_t vq = vp.unaryExpr([](auto x){ return std::max(x, 0.); }); 
    vq *= ni / vq.sum();

    vec_t ww(no); 
    bvec_t ju(ni, false);
    vec_t xm(ni); 
    vec_t xv;
    vec_t xs(ni); xs.setZero();
    if (kopt == 2) { xv.setZero(ni); }

    SpChkvars::eval(x,ju);
    if (jd(0) > 0) {
        for (int i = 1; i < jd(0) + 1; ++i) {
            ju[jd(i)-1] = false;
        }
    }

    // TODO: this error should be generalized. Gaussian too.
    // can't find true value in ju
    if (std::find_if(ju.begin(), ju.end(), [](auto x) { return x;}) == ju.end()) {
        jerr=7777; 
        return;
    } 

    for (int i = 0; i < no; ++i) {
        ww(i) = y.row(i).sum();
        if (ww(i)) y.row(i) /= ww(i);
    }
    auto sw = ww.sum();
    ww /= sw;

    if (nc != 1 && kopt == 2) { 
        MultSpLStandardize2::eval(x,ww,ju,isd,intr,xm,xs,xv);
    }
    else {
        SpLStandardize2::eval(x,ww,ju,isd,intr,xm,xs);
    }

    if (isd) { 
        for (int j = 0; j < ni; ++j) {
            cl.col(j) *= xs(j);
        }
    }

    if (nc == 1) { 
        auto y_1 = y.col(0);
        auto g_1 = g.col(0);
        Eigen::Map<Eigen::MatrixXd> ca_slice(ca.data(), nx, nlam);
        Eigen::Map<Eigen::VectorXd> a0_slice(a0.data(), a0.size());
        splognet2n<FloatType, apply_suggested_fixes>(
                parm,x,y_1,g_1,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
             thr,isd,intr,maxit,kopt,xm,xs,lmu,a0_slice,ca_slice,ia,nin,dev0,dev,alm,nlp,jerr);
    }
    else if (kopt == 2) {
        multsplognetn<FloatType, apply_suggested_fixes>(
                parm,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
                thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
    }
    else { 
        splognetn<FloatType, apply_suggested_fixes>(
                parm,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,thr,
             isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
    }

    if(jerr > 0) return; 
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
            if (!intr) { a0(ic,k)=0.0; }
            else { 
                for (int i = 0; i < nk; ++i) {
                    a0(ic, k) -= ca_slice(i, ic) * xm(ia(i)-1);
                }
            }
        }
    }
}

} // namespace transl
} // namespace glmnetpp

