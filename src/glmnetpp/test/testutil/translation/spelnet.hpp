#pragma once
#include <Eigen/Core>
#include <vector>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <testutil/translation/spelnet1.hpp>
#include <testutil/translation/spelnet2.hpp>

namespace glmnetpp {
namespace transl {

template <class FloatType
        , class ValueType
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
        , class AlmType>
inline void spelnet(
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
        IntType& jerr
            )
{
    using value_t = ValueType;
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using bvec_t = std::vector<bool>;

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

    SpChkvars::eval(x, ju);

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

    if (!ka) {
        g.setZero(ni);
        SpStandardize::eval(x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv);
    }
    else {
        SpStandardize1::eval(x,y,w,isd,intr,ju,xm,xs,ym,ys,xv);
    }

    cl /= ys; 
    if (isd) { 
        for (int j = 0; j < ni; ++j) {
            cl.col(j) *= xs(j);
        }
    }

    if (flmin >= 1.0) vlam = ulam / ys;

    // ka == 0
    if (!ka) {
        spelnet1<FloatType>(
                parm,ju,vp,cl,g,w,ne,nx,x,nlam,flmin,vlam,thr,maxit,
                xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
    }
    else {
        spelnet2<FloatType>(
                parm,ju,vp,cl,y,w,ne,nx,x,nlam,flmin,vlam,thr,maxit,
                xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
    }

    if (jerr > 0) return;

    for (int k = 0; k < lmu; ++k) {
        alm(k) *= ys; 
        auto nk = nin(k);
        for (int l = 0; l < nk; ++l) {
            ca(l,k) *= ys / xs(ia(l));
        }
        a0(k)=0.0;
        if (intr) {
            for (int i = 0; i < nk; ++i) {
                a0(k) -= ca(i, k) * xm(ia(i));
            }
            a0(k) += ym;
        }
    }
}

} // namespace transl
} // namespace glmnetpp
