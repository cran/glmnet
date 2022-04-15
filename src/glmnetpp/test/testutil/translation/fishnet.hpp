#pragma once
#include <Eigen/Core>
#include <vector>

// have been separately unittested
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>        
#include <glmnetpp_bits/elnet_driver/standardize.hpp> 
#include <testutil/translation/fishnet1.hpp>

namespace glmnetpp {
namespace transl {

template <class FloatType
        , class ValueType
        , class XType
        , class YType
        , class GType
        , class WType
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
inline void fishnet(
        ValueType parm,
        XType& x,
        YType& y,
        GType& g,
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
        bool isd,
        bool intr,
        IntType maxit,
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
    using bvec_t = std::vector<bool>;

    if (vp.maxCoeff() <= 0) { jerr = 10000; return; }
    if (y.minCoeff() < 0) { jerr = 8888; return; }

    auto no = x.rows();
    auto ni = x.cols();
    vec_t vq = vp.unaryExpr([](auto x){ return std::max(x, 0.); }); 
    vq *= ni / vq.sum();

    vec_t ww(no); 
    bvec_t ju(ni, false);
    vec_t xm(ni); 
    vec_t xs;
    if (isd) { xs.setZero(ni); }

    Chkvars::eval(x,ju);
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

    ww = w.array().max(0.0);
    auto sw = ww.sum();
    if (sw <= 0.0) { jerr = 9999; return; }
    ww /= sw;

    LStandardize1::eval(x,ww,ju,isd,intr,xm,xs);
    if (isd) { 
        for (int j = 0; j < ni; ++j) {
            cl.col(j) *= xs(j);
        }
    }

    fishnet1<FloatType>(
            parm,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
            thr,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);

    if(jerr > 0) return; 
    dev0 *= 2.0*sw;
    for (int k = 0; k < lmu; ++k) {
        auto nk = nin(k);
        if (isd) { 
            for (int l = 0; l < nk; ++l) {
                ca(l, k) /= xs(ia(l)-1);
            }
        }
        if (!intr) { a0(k)=0.0; }
        else { 
            for (int i = 0; i < nk; ++i) {
                a0(k) -= ca(i, k) * xm(ia(i)-1);
            }
        }
    }
}

} // namespace transl
} // namespace glmnetpp

