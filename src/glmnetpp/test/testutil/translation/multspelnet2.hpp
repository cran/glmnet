#pragma once
#include <cmath>
#include <testutil/translation/multelnet2.hpp>
#include <Eigen/Core>

namespace glmnetpp {
namespace transl {

template <class FloatType
        , class IntType
        , class AType
        , class IAType
        , class XType
        , class YType
        , class WType
        , class XMType
        , class XSType
        , class XVType
        , class VPType
        , class ValueType
        , class CLType
        , class GJType
        , class GKType
        , class ISCType
        , class DelType
        , class OType>
inline void multspelnet2_do_b(
    bool& iz,
    bool& jz,
    IntType nin,
    AType& a,
    const IAType& ia,
    IntType& nlp,
    const XType& x,
    YType& y,
    const WType& w,
    const XMType& xm,
    const XSType& xs,
    const XVType& xv,
    const VPType& vp,
    ValueType ab,
    ValueType dem,
    const CLType& cl,
    GJType& gj,
    GKType& gk,
    ISCType& isc,
    DelType& del,
    OType& o,
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
            auto ak = a.col(k);
            auto nr = ak.size();

            // TODO
            for (int j = 0; j < nr; ++j) {
                gj(j) = x.col(k).cwiseProduct(w).dot(
                        (y.col(j).array() + o(j)).matrix()) / xs(k);
            }

            gk = gj + xv(k) * a.col(k);
            auto gkn = gk.norm();
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
            del.array() = ak.array() - del.array();
            if (del.array().abs().maxCoeff() <= 0.0) continue;

            rsq -= (del.array() * (2.0 * gj.array() - xv(k) * del.array())).sum();
            for (int j = 0; j < y.cols(); ++j) { y.col(j) -= (del(j)/xs(k)) * x.col(k); }
            // TODO
            o += (xm(k)/xs(k)) * del;
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
        , class WType
        , class IntType
        , class XType
        , class ULamType
        , class XMType
        , class XSType
        , class XVType
        , class AOType
        , class IAType
        , class KinType
        , class RsqoType
        , class AlmoType
        >
inline void multspelnet2(
        ValueType beta,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        YType& y,
        const WType& w,
        IntType ne,
        IntType nx,
        const XType& x,
        IntType nlam,
        ValueType flmin,
        const ULamType& ulam,
        ValueType thr,
        IntType maxit,
        const XMType& xm,
        const XSType& xs,
        const XVType& xv,
        ValueType& ys0,
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
    Eigen::VectorXd o(nr); o.setZero(); // TODO: mean shifts

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
        // TODO: make faster
        // also, not sure if I can do this...
        gj = y.transpose() * x.col(j).cwiseProduct(w);
        g(j) = gj.norm() / xs(j);
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
            if (iz * jz != 0) { 
                multspelnet2_do_b<FloatType>(
                        iz, jz, nin, a, ia, nlp, x, y, w, xm, xs, xv, vp,
                        ab, dem, cl, gj, gk, isc, del, o, rsq, thr, maxit);
            }
            while (1) {

                bool converged_kkt = false;
                // :again:
                while (1) {
                    if (nlp > maxit) throw util::maxit_reached_error();
                    ++nlp;
                    auto dlx = 0.0;

                    for (int k = 0; k < ni; ++k) {
                        if (ix(k) == 0) continue;
                        // TODO:
                        for (int j = 0; j < nr; ++j) {
                            gj(j) = x.col(k).cwiseProduct(w).dot(
                                    (y.col(j).array() + o(j)).matrix()) / xs(k);
                        }
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
                        for (int j = 0; j < nr; ++j) { y.col(j) -= (del(j) / xs(k)) * x.col(k); }
                        // TODO:
                        o += (xm(k) / xs(k)) * del;
                        dlx = std::max(dlx, xv(k) * del.array().square().maxCoeff());
                    }

                    if (dlx < thr) {
                        bool ixx = false;
                        for (int k = 0; k < ni; ++k) {
                            if (ix(k) == 1) continue;
                            if (ju[k] == 0) continue;
                            // TODO: since it's matrix product,
                            // it might create temporary
                            for (int j = 0; j < nr; ++j) {
                                gj(j) = x.col(k).cwiseProduct(w).dot(
                                        (y.col(j).array() + o(j)).matrix()) / xs(k);
                            }
                            g(k) = gj.norm();
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

                multspelnet2_do_b<FloatType>(
                        iz, jz, nin, a, ia, nlp, x, y, w, xm, xs, xv, vp,
                        ab, dem, cl, gj, gk, isc, del, o, rsq, thr, maxit);
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
