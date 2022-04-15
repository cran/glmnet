#pragma once
#include <Eigen/Core>
#include <testutil/mock_pb.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>

namespace glmnetpp {
namespace transl {

// Note: the FloatType is really important!
// In Fortran, constants are by default single-precision,
// so to test compatability, we have set the FloatType = float.
// But to test against our abstracted version of elnet1,
// we should use FloatType = double.

// Brute-force translation of legacy elnet1 Fortran.
template <class FloatType
        , class IntType
        , class DAType
        , class AType
        , class IAType
        , class GType
        , class XVType
        , class VPType
        , class ValueType
        , class JUType
        , class CLType
        , class CType
        , class MMType>
inline void spelnet1_do_b(
    bool& iz,
    bool& jz,
    IntType ni,
    IntType nin,
    IntType m,
    DAType& da,
    AType& a,
    const IAType& ia,
    IntType& nlp,
    GType& g,
    const XVType& xv,
    const VPType& vp,
    ValueType ab,
    ValueType dem,
    const JUType& ju,
    const CLType& cl,
    ValueType& rsq,
    const CType& c,
    const MMType& mm,
    ValueType thr,
    IntType maxit,
    IntType& jerr
        )
{
    using int_t = IntType;

    iz = true;
    for (int_t j = 0; j < nin; ++j) {
        da(j) = a(ia(j));
    }
    while (1) {
        ++nlp;
        auto dlx = 0.0;
        for (int_t l = 0; l < nin; ++l) {
            auto k = ia(l);
            auto ak = a(k);
            auto u = g(k) + ak * xv(k);
            auto v = std::abs(u) - vp(k) * ab;
            a(k) = 0.0;
            if (v > 0.0) {
                a(k) = std::max(cl(0, k), 
                        std::min(cl(1, k), 
                            std::copysign(v,u)/(xv(k)+vp(k)*dem)) );
            }
            if (a(k) == ak) continue;

            auto del = a(k) - ak;
            rsq += del * (static_cast<FloatType>(2.0) * g(k) - del * xv(k));
            dlx = std::max(xv(k) * del * del, dlx);
            for (int_t j = 0; j < nin; ++j) {
                g(ia(j)) -= c(ia(j), mm(k)-1) * del;
            }
        }
        if (dlx < thr) break;
        if (nlp > maxit) { 
            jerr = -m-1; 
            throw std::exception();
        }
    }

    for (int_t j = 0; j < nin; ++j) {
        da(j) = a(ia(j)) - da(j);
    }

    for (int_t j = 0; j < ni; ++j) {
        if (mm(j) != 0) continue;
        if (ju[j] != 0) {
            g(j) -= da.head(nin).dot(c.row(j).head(nin));
        }
    }
    jz = false;
}

template <class FloatType
        , class ValueType
        , class JUType
        , class VPType
        , class CLType
        , class GType
        , class WType
        , class IntType
        , class XMType
        , class XSType
        , class XType
        , class ULamType
        , class XVType
        , class AOType
        , class IAType
        , class KinType
        , class RSQOType
        , class ALMOType>
inline void spelnet1(
    ValueType beta,
    const JUType& ju,
    const VPType& vp,
    const CLType& cl,
    GType& g,
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

    Eigen::MatrixXd c(ni, nx); c.setZero();
    Eigen::VectorXd a(ni); a.setZero();
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXd da(ni); da.setZero();

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

    for (int_t m = 0; m < nlam; ++m) {

        if (int_param_t::itrace != 0) mock_setpb(m);
        if (flmin >= 1.0) { alm = ulam(m); }
        else if (m > 1) { alm *= alf; }
        else if (m == 0) { alm = int_param_t::big; }
        else {
            alm = 0.0;
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j] == 0) continue;
                if (vp(j) <= 0.0) continue;
                alm = std::max(alm, std::abs(g(j)) / vp(j));
            }
            alm = alf * alm / std::max(bta, 1e-3);
        }
        auto dem = alm * omb;
        auto ab = alm * bta;
        auto rsq0 = rsq;

        try {

            bool jz = true;

            while (1) {
                if (iz * jz != 0) { 
                    spelnet1_do_b<FloatType>(
                            iz, jz, ni, nin, m, da, a, ia, nlp, g,
                            xv, vp, ab, dem, ju, cl, rsq, c, mm, thr, maxit, jerr); 
                }

                ++nlp;
                auto dlx = 0.0;
                for (int_t k = 0; k < ni; ++k) {
                    if (ju[k] == 0) continue;
                    auto ak = a(k);
                    auto u = g(k) + ak * xv(k);
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
                        for (int_t j = 0; j < ni; ++j) {
                            if (ju[j] == 0) continue;
                            if (mm(j) != 0) { 
                                c(j, nin-1) = c(k, mm(j)-1);
                                continue;
                            }
                            if (j == k) {
                                c(j, nin-1) = xv(j);
                                continue;
                            }
                            auto wx_k = x.col(k).cwiseProduct(w);
                            c(j,nin-1) = 
                                 (x.col(j).dot(wx_k) - xm(j) * xm(k))
                                 / (xs(j) * xs(k));
                        }
                        mm(k) = nin;
                        ia(nin-1) = k;
                    }

                    auto del = a(k) - ak;
                    rsq += del * (static_cast<FloatType>(2.0) * g(k) - del * xv(k));
                    dlx = std::max(xv(k) * del * del, dlx);
                    for (int_t j = 0; j < ni; ++j) {
                        if (ju[j] != 0) {
                            g(j) -= c(j, mm(k)-1) * del;
                        }
                    }
                }

                if (dlx < thr) break;
                if (nin > nx) break;
                if (nlp > maxit) { 
                    jerr = -m-1; 
                    return;
                }

                spelnet1_do_b<FloatType>(
                        iz, jz, ni, nin, m, da, a, ia, nlp, g,
                        xv, vp, ab, dem, ju, cl, rsq, c, mm, thr, maxit, jerr); 
            }

        }
        catch (const std::exception&) {
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
            if (ao(j,m) != 0.0) ++me;
        }
        if (me > ne) break;
        if (rsq-rsq0 < int_param_t::sml*rsq) break;
        if (rsq > int_param_t::rsqmax) break;
    }
}

} // namespace transl
} // namespace glmnetpp
