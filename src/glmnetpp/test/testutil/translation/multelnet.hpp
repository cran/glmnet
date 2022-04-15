#pragma once
#include <vector>
#include <testutil/translation/multelnet2.hpp>
#include <Eigen/Core>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

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
        , class IntType
        , class ULamType
        , class A0Type
        , class AOType
        , class IAType
        , class KinType
        , class RsqoType
        , class AlmoType
        >
inline void multelnet(
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
        IntType& jerr
        )
{
    if (vp.maxCoeff() <= 0.0) { jerr=10000; return; }
    auto ni = x.cols();
    auto nr = y.cols();
    Eigen::VectorXd vq=vp.array().max(0.0).matrix(); 
    vq*=ni/vq.sum();

    Eigen::VectorXd clt(2 * nr * ni);
    Eigen::VectorXd xm(ni);
    Eigen::VectorXd xs(ni);
    Eigen::VectorXd ym(nr);
    Eigen::VectorXd ys(nr);
    std::vector<bool> ju(ni, false);
    Eigen::VectorXd xv(ni);
    
    Chkvars::eval(x, ju);

    if (jd(0) > 0) {
        for (int i = 1; i < jd(0) + 1; ++i) {
            ju[jd(i)-1] = false;
        }
    }
    if (std::find_if(ju.begin(), ju.end(), [](auto x) { return x;}) == ju.end()) {
        jerr=7777; 
        return;
    } 

    double ys0 = 0.0;
    MultStandardize1::eval(x, y, w, isd, jsd, intr, ju, xm, xs, ym, ys, xv, ys0);

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

    multelnet2<FloatType>(
            beta, ju, vq, clt, y, ne, nx, x, nlam, flmin, ulam, thr, maxit, xv,
            ys0, lmu, ca, ia, nin, rsq, alm, nlp, jerr);

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

} // namespace transl
} // namespace glmnetpp
