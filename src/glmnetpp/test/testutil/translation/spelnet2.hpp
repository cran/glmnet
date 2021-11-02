#pragma once
#include <Eigen/Core>
#include <testutil/mock_pb.hpp>
#include <glmnetpp_bits/internal.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>

namespace glmnetpp {
namespace transl {

// Note: the FloatType is really important!
// In Fortran, constants are by default single-precision,
// so to test compatability, we have set the FloatType = float.
// But to test against our abstracted version of elnet1,
// we should use FloatType = double.

// Brute-force translation of legacy elnet2 Fortran.

template <class FloatType
        , class IntType
        , class AType
        , class IAType
        , class WType
        , class XType
        , class YType
        , class XMType
        , class XSType
        , class XVType
        , class VPType
        , class ValueType
        , class CLType>
inline void spelnet2_do_b(
    bool& iz,
    bool& jz,
    IntType nin,
    IntType m,
    AType& a,
    const IAType& ia,
    const WType& w,
    IntType& nlp,
    const XType& x,
    YType& y,
    const XMType& xm,
    const XSType& xs,
    const XVType& xv,
    const VPType& vp,
    ValueType ab,
    ValueType dem,
    const CLType& cl,
    ValueType& o,
    ValueType& rsq,
    ValueType thr,
    IntType maxit,
    IntType& jerr
        )
{
    using int_t = IntType;

    iz = true;
    while (1) {
        ++nlp;
        auto dlx = 0.0;
        for (int_t l = 0; l < nin; ++l) {
            auto k = ia(l);
            auto gk = x.col(k).cwiseProduct(w).dot(
                        (y.array() + o).matrix()) / xs(k);
            auto ak = a(k);
            auto u = gk + ak * xv(k);
            auto v = std::abs(u) - vp(k) * ab;
            a(k) = 0.0;
            if (v > 0.0) {
                a(k) = std::max(cl(0, k), 
                        std::min(cl(1, k), 
                            std::copysign(v,u)/(xv(k)+vp(k)*dem)) );
            }
            if (a(k) == ak) continue;

            auto del = a(k) - ak;
            rsq += del * (static_cast<FloatType>(2.0) * gk - del * xv(k));
            auto del_scaled = del / xs(k);
            y -= del_scaled * x.col(k); 
            o += del_scaled * xm(k);
            dlx = std::max(xv(k) * del * del, dlx);
        }
        if (dlx < thr) break;
        if (nlp > maxit) { 
            jerr = -m-1; 
            throw util::maxit_reached_error();
        }
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
        , class RSQOType
        , class ALMOType>
inline void spelnet2(
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
    IntType& lmu,
    AOType& ao,
    IAType& ia,
    KinType& kin,
    RSQOType& rsqo,
    ALMOType& almo,
    IntType& nlp,
    IntType& jerr) 
{
    using int_t = IntType;
    using int_param_t = InternalParams;

    int_t ni = x.cols();

    Eigen::VectorXd a(ni); a.setZero();
    Eigen::VectorXd g(ni); g.setZero();
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXi ix(ni); ix.setZero();

    auto bta = beta;
    auto omb = 1.0 - beta;
    auto alm = 0.0;
    auto alf = 1.0;

    if (flmin < 1.0) {
        auto eqs = std::max(int_param_t::eps, flmin);
        alf = std::pow(eqs, static_cast<FloatType>(1.0)/(nlam - 1));
    }

    auto rsq = 0.0;
    nlp = 0;
    int_t nin = 0;
    bool iz = false;
    auto mnl = std::min(int_param_t::mnlam, nlam);

    auto o = 0.0;

    for (int_t j = 0; j < ni; ++j) {
        if (ju[j] == 0) continue;
        g(j) = std::abs(
                x.col(j).cwiseProduct(w).dot(
                    (y.array() + o).matrix()) / xs(j));
    }

    for (int_t m = 0; m < nlam; ++m) {

        if (int_param_t::itrace != 0) mock_setpb(m);
        auto alm0 = alm;
        if (flmin >= 1.0) { alm = ulam(m); }
        else if (m > 1) { alm *= alf; }
        else if (m == 0) { alm = int_param_t::big; }
        else {
            alm0 = 0.0;
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j] == 0) continue;
                if (vp(j) > 0.0) {
                    alm0 = std::max(alm0, g(j) / vp(j));
                }
            }
            alm0 /= std::max(bta, 1e-3);
            alm = alf * alm0;
        }
        auto dem = alm * omb;
        auto ab = alm * bta;
        auto rsq0 = rsq;

        try {

            bool jz = true;
            auto tlam = bta * (static_cast<FloatType>(2.0) * alm - alm0);
            for (int_t k = 0; k < ni; ++k) {
                if (ix(k) == 1) continue;
                if (ju[k] == 0) continue;
                if (g(k) > tlam * vp(k)) ix(k) = 1;
            }

            while (1) {

                if (iz * jz != 0) { 
                    spelnet2_do_b<FloatType>(
                            iz, jz, nin, m, a, ia, w, nlp, x, y,
                            xm, xs, xv, vp, ab, dem, cl, o, rsq, thr, maxit, jerr); 
                }


                // :again:
                bool converged_kkt = false;
                while (1) {

                    if (nlp > maxit) {
                        jerr = -m-1;
                        return;
                    }

                    ++nlp;
                    auto dlx = 0.0;
                    for (int_t k = 0; k < ni; ++k) {
                        if (ix(k) == 0) continue;
                        auto gk = x.col(k).cwiseProduct(w).dot(
                                    (y.array() + o).matrix()) / xs(k);
                        auto ak = a(k);
                        auto u = gk + ak * xv(k);
                        auto v = std::abs(u) - vp(k) * ab;
                        a(k) = 0.0;
                        if (v > 0.0) {
                            a(k) = std::max(cl(0, k), 
                                    std::min(cl(1, k), 
                                        std::copysign(v,u)/(xv(k)+vp(k)*dem)) );
                        }
                        if (a(k) == ak) continue;
                        if (mm(k) == 0) {
                            ++nin;
                            if (nin > nx) break;
                            mm(k) = nin;
                            ia(nin-1) = k;
                        }

                        auto del = a(k) - ak;
                        rsq += del * (static_cast<FloatType>(2.0) * gk - del * xv(k));
                        auto del_scaled = del / xs(k);
                        y -= del_scaled * x.col(k); 
                        o += del_scaled * xm(k);
                        dlx = std::max(xv(k) * del * del, dlx);
                    }

                    if (nin > nx) throw util::max_active_reached_error();
                    if (dlx < thr) {
                        bool ixx = false;
                        for (int_t k = 0; k < ni; ++k) {
                            if (ix(k) == 1) continue;
                            if (ju[k] == 0) continue;
                            g(k) = std::abs(
                                    x.col(k).cwiseProduct(w).dot(
                                        (y.array() + o).matrix()) / xs(k));
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

                if (nlp > maxit) { 
                    jerr = -m-1; 
                    return;
                }

                spelnet2_do_b<FloatType>(
                        iz, jz, nin, m, a, ia, w, nlp, x, y,
                        xm, xs, xv, vp, ab, dem, cl, o, rsq, thr, maxit, jerr); 
            }

        }
        catch (const util::max_active_reached_error&) {}
        catch (const util::maxit_reached_error&) {
            return;
        }

        if (nin > nx) { jerr = -10000-m-1; break; }
        if (nin > 0) { 
            for (int_t j = 0; j < nin; ++j) {
                ao(j, m) = a(ia(j));
            }
        }
        kin(m) = nin;
        rsqo(m) = rsq;
        almo(m) = alm;
        lmu = m+1;
        if (lmu < mnl) continue;
        if (flmin >= 1.0) continue;
        auto me = 0;
        for (int_t j = 0; j < nin; ++j) {
            if (ao(j,m)) ++me;
        }
        if (me > ne) break;
        if (rsq-rsq0 < int_param_t::sml*rsq) break;
        if (rsq > int_param_t::rsqmax) break;
    }
}

} // namespace transl
} // namespace glmnetpp
