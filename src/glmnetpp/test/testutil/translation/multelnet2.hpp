#pragma once
#include <cmath>
#include <testutil/mock_pb.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <Eigen/Core>

namespace glmnetpp {
namespace transl {

inline double bnorm(
        double b0,
        double al1p,
        double al2p,
        double g,
        double usq
        )
{
    constexpr double thr = 1e-10;
    constexpr int mxit = 100;
    auto b = b0;
    auto zsq = b*b + usq;
    if (zsq <= 0.0) { return 0.0; }
    auto z = std::sqrt(zsq);
    auto f = b*(al1p+al2p/z)-g;
    int it = 0;
    for (; it < mxit; ++it) {
        b -= f/(al1p+al2p*usq/(z*zsq));
        zsq = b*b + usq;
        if (zsq <= 0.0) { return 0.0; }
        z = std::sqrt(zsq);
        f = b*(al1p+al2p/z)-g;
        if (std::abs(f) <= thr) break;
        if (b <= 0.0) { b = 0.0; break; }
    }
    if (it >= mxit) throw util::bnorm_maxit_reached_error();
    return b;
}

template <class GKType
        , class ValueType
        , class CLType
        , class AType
        , class ISCType>
inline void chkbnds(
        const GKType& gk,
        ValueType gkn,
        ValueType xv,
        const CLType& cl,
        ValueType al1,
        ValueType al2,
        AType& a,
        ISCType& isc
        ) 
{
    auto al1p = 1.0 + al1/xv;
    auto al2p = al2/xv;
    isc.setZero();
    auto gsq = gkn*gkn;
    auto asq = a.squaredNorm();
    double usq = 0.0;
    double u = 0.0;
    int kn = -1;
    while (1) {
        double vmx = 0.0;
        for (int k = 0; k < a.size(); ++k) {
            auto v = std::max(a(k)-cl(1,k), cl(0,k)-a(k));
            if (v > vmx) { vmx = v; kn = k; }
        }
        if (vmx <= 0.0) break;
        if (isc(kn)) break;
        gsq -= gk(kn)*gk(kn);
        auto g = std::sqrt(std::abs(gsq))/xv;
        if (a(kn) < cl(0, kn)) u = cl(0, kn);
        if (a(kn) > cl(1, kn)) u = cl(1, kn);
        usq += u*u;
        double b = 0.0;
        if (usq == 0.0) b = std::max(0., (g-al2p)/al1p);
        else {
            auto b0 = std::sqrt(asq - a(kn) * a(kn));
            b = bnorm(b0, al1p, al2p, g, usq);
        }
        asq = usq + b*b;
        if (asq <= 0.0) { a.setZero(); break; }
        a(kn) = u;
        isc(kn) = 1;
        auto f = 1.0/(xv * (al1p+al2p/std::sqrt(asq)));
        for (int j = 0; j < a.size(); ++j) {
            if (!isc(j)) a(j) = f * gk(j);
        }
    }
}

template <class FloatType
        , class IntType
        , class AType
        , class IAType
        , class XType
        , class YType
        , class XVType
        , class VPType
        , class ValueType
        , class CLType
        , class GJType
        , class GKType
        , class ISCType
        , class DelType>
inline void multelnet2_do_b(
    bool& iz,
    bool& jz,
    IntType nin,
    AType& a,
    const IAType& ia,
    IntType& nlp,
    const XType& x,
    YType& y,
    const XVType& xv,
    const VPType& vp,
    ValueType ab,
    ValueType dem,
    const CLType& cl,
    GJType& gj,
    GKType& gk,
    ISCType& isc,
    DelType& del,
    ValueType& rsq,
    ValueType thr,
    IntType maxit
        )
{
    iz = true;
    while (1) {
        ++nlp;
        auto dlx = 0.0;
        for (int l = 0; l < nin; ++l) {
            auto k = ia(l)-1;
            // TODO
            gj = y.transpose() * x.col(k);
            gk = gj + xv(k) * a.col(k);
            auto gkn = gk.norm();
            auto ak = a.col(k);
            auto u = 1.0 - ab*vp(k)/gkn;
            del = ak;
            if (u <= 0.0) { ak.setZero(); }
            else {
                ak = gk*(u/(xv(k)+dem*vp(k)));

                auto nr = ak.size();
                Eigen::Map<const Eigen::MatrixXd> cl_slice(
                        cl.data() + k * 2 * nr, 2, nr);

                chkbnds(gk,gkn,xv(k),cl_slice,
                        dem*vp(k),ab*vp(k),ak,isc);
            }
            del.array() = ak.array() - del.array();
            if (del.array().abs().maxCoeff() <= 0.0) continue;

            rsq -= (del.array() * (2.0 * gj.array() - xv(k) * del.array())).sum();
            for (int j = 0; j < y.cols(); ++j) { y.col(j) -= del(j) * x.col(k); }
            dlx = std::max(dlx, xv(k) * del.array().square().maxCoeff());
        }
        if (dlx < thr) break;
        if (nlp > maxit) throw util::maxit_reached_error();
    }

    jz = false;
}

template <class FloatType
        , class ValueType
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
        , class RsqoType
        , class AlmoType
        >
inline void multelnet2(
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
        RsqoType& rsqo,
        AlmoType& almo,
        IntType& nlp,
        IntType& jerr
        )
{
    using int_param_t = InternalParams;

    auto nr = y.cols();
    auto ni = x.cols();
    auto eps = int_param_t::eps;
    auto sml = int_param_t::sml;
    auto rsqmax = int_param_t::rsqmax;
    auto mnlam = int_param_t::mnlam;
    auto itrace = int_param_t::itrace;
    auto big = int_param_t::big;

    Eigen::MatrixXd a(nr, ni); a.setZero();
    Eigen::VectorXd gj(nr);
    Eigen::VectorXd gk(nr);
    Eigen::VectorXd del(nr);
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXd g(ni);
    Eigen::VectorXi ix(ni); ix.setZero(); // strong set?
    Eigen::VectorXi isc(nr);

    auto bta = beta;
    auto omb = 1.0 - bta;

    thr *= ys0 / nr;
    auto alf = 1.0;

    if (flmin < 1.0) { 
        auto eqs=std::max(eps,flmin); 
        alf=std::pow(eqs, (1.0/(nlam-1))); 
    }
    auto rsq = ys0; 
    nlp = 0;
    auto nin =0; 
    bool iz=false; 
    auto mnl=std::min(mnlam,nlam); 
    auto alm=0.0;
    for (int j = 0; j < ni; ++j) {
        if (!ju[j]) continue;
        // TODO
        g(j) = (y.transpose() * x.col(j)).norm();
    }

    for (int m = 0; m < nlam; ++m) {
        if (itrace) mock_setpb(m); 
        auto alm0=alm;
        if (flmin >= 1.0) { alm=ulam(m); }
        else if (m > 1) { alm*=alf; }
        else if (m == 0) { alm=big; }
        else { 
            alm0=0.0;
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue;
                if (vp(j) > 0.0) alm0=std::max(alm0,g(j)/vp(j));
            }
            alm0 /= std::max(bta,1e-3); 
            alm=alf*alm0;
        }
        auto dem=alm*omb; 
        auto ab=alm*bta; 
        auto rsq0=rsq; 
        Eigen::Map<Eigen::MatrixXd> ao_slice(
                ao.data() + m * nx * nr, nx, nr);
        bool jz=1;
        auto tlam=bta*(2.0*alm-alm0);
        for (int k = 0; k < ni; ++k) {
            if (ix(k) || !ju[k]) continue; 
            if (g(k) > tlam*vp(k)) ix(k)=1;
        }

        try {
            while (1) {
                if (iz * jz != 0) { 
                    multelnet2_do_b<FloatType>(
                            iz, jz, nin, a, ia, nlp, x, y, xv, vp,
                            ab, dem, cl, gj, gk, isc, del, rsq, thr, maxit);
                }

                bool converged_kkt = false;
                // :again:
                while (1) {
                    if (nlp > maxit) throw util::maxit_reached_error();
                    ++nlp;
                    auto dlx = 0.0;

                    for (int k = 0; k < ni; ++k) {
                        if (ix(k) == 0) continue;
                        gj = y.transpose() * x.col(k);
                        gk = gj + xv(k) * a.col(k);
                        auto gkn = gk.norm();
                        auto ak = a.col(k);
                        auto u = 1.0 - ab*vp(k)/gkn;
                        del = ak;
                        if (u <= 0.0) { ak.setZero(); }
                        else {
                            ak = gk*(u/(xv(k)+dem*vp(k)));
                            Eigen::Map<const Eigen::MatrixXd> cl_slice(
                                    cl.data() + k * 2 * nr, 2, nr);
                            chkbnds(gk,gkn,xv(k),cl_slice,
                                    dem*vp(k),ab*vp(k),ak,isc);
                        }
                        del = ak - del;
                        if (del.array().abs().maxCoeff() <= 0.0) continue;
                        if (mm(k) == 0) {
                            ++nin;
                            if (nin > nx) throw util::max_active_reached_error();
                            mm(k) = nin;
                            ia(nin-1) = k+1;
                        }

                        rsq -= (del.array() * (2.0 * gj - xv(k) * del).array()).sum();
                        for (int j = 0; j < nr; ++j) { y.col(j) -= del(j) * x.col(k); }
                        dlx = std::max(dlx, xv(k) * del.array().square().maxCoeff());
                    }

                    if (dlx < thr) {
                        bool ixx = false;
                        for (int k = 0; k < ni; ++k) {
                            if (ix(k) == 1) continue;
                            if (ju[k] == 0) continue;
                            // TODO: since it's matrix product,
                            // it might create temporary
                            g(k) = (y.transpose() * x.col(k)).norm();
                            if (g(k) > ab * vp(k)) {
                                ix(k) = 1;
                                ixx = true;
                            }
                        }
                        if (ixx) continue;
                        converged_kkt = true; 
                    }
                    break;
                }

                if (converged_kkt) break; 
                if (nlp > maxit) throw util::maxit_reached_error();

                multelnet2_do_b<FloatType>(
                        iz, jz, nin, a, ia, nlp, x, y, xv, vp,
                        ab, dem, cl, gj, gk, isc, del, rsq, thr, maxit);
            }
        } 
        catch (const util::maxit_reached_error& e) {
            jerr = e.err_code(m);
            return;
        }
        catch (const util::bnorm_maxit_reached_error& e) {
            jerr = e.err_code(m);
            return;
        }
        catch (const util::max_active_reached_error& e) {
            jerr = e.err_code(m);
            break;
        }

        for (int j = 0; j < nr; ++j) {
            for (int i = 0; i < nin; ++i) {
                ao_slice(i,j) = a(j, ia(i)-1);
            }
        }
        kin(m)=nin;
        rsqo(m)=1.0-rsq/ys0; 
        almo(m)=alm; 
        lmu=m+1;
        if(lmu < mnl || flmin >= 1.0) continue;
        auto me = (ao_slice.col(0).head(nin).array() != 0.0).count();
        if(me > ne ||
           rsq0 - rsq < sml * rsq ||
           rsqo(m) > rsqmax) break;
    }
}

} // namespace transl
} // namespace glmnetpp
