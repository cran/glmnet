#pragma once
#include <Eigen/Core>
#include <testutil/mock_pb.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>

namespace glmnetpp {
namespace transl {

template <class YType
        , class GType
        , class QType
        , class AZType>
inline auto kazero(
        const YType& y,
        const GType& g,
        const QType& q,
        AZType&& az
        )
{
    auto n = y.rows();
    auto kk = y.cols();
    Eigen::MatrixXd e(n, kk);
    Eigen::VectorXd s(n);
    az.setZero();
    e.array() = g.array().exp();
    s = e.rowwise().sum();
    double dm;
    do {
        dm = 0.0;
        for (int k = 0; k < kk; ++k) {
            auto t = 0.0;
            auto u = 0.0;
            for (int i = 0; i < n; ++i) {
                auto pik = e(i,k)/s(i);
                t += q(i) * (y(i,k) - pik);
                u += q(i) * pik * (1.0-pik);
            }
            auto d = t/u;
            az(k) += d;
            auto ed = std::exp(d);
            dm = std::max(dm, std::abs(d));
            for (int i = 0; i < n; ++i) {
                auto z = e(i,k);
                e(i,k) = z * ed;
                s(i) += -z + e(i,k);
            }
        }
    } while(dm >= 1e-7);
    az.array() -= az.sum() / kk;
}

template <class BetaType
        , class CLType
        , class AType
        , class MType>
inline auto elc(
        BetaType parm,
        const CLType& cl,
        const AType& a,
        MType&& m
        )
{
    auto n = a.size();
    auto fn = n; 
    auto am = a.sum()/fn;
    auto out = am;
    if (parm && (n != 2)) { 
        // Note: it is VERY important that we take head.
        // Otherwise, m is resized, which will make it undefined behavior for the caller.
        // We assume m has already been reshaped such that its size is at least n.
        m.head(n) = Eigen::VectorXi::LinSpaced(n, 0, n-1);
        std::sort(m.data(), m.data() + n, 
                [&](int i, int j) { return a[i] < a[j]; });

        if (a(m(0)) == a(m(n-1))) {
            out = a(0); 
        } else {
            double ad = 0.0;
            if (n % 2 == 1) { ad = a(m(n/2)); }
            else { ad = 0.5*(a(m(n/2))+a(m(n/2-1))); }
            if (parm == 1.0) { 
                out = ad; 
            } else {
                auto b1 = std::min(am, ad); 
                auto b2 = std::max(am, ad); 
                auto k2=1;
                while (a(m(k2-1)) <= b1) { ++k2; }
                auto k1 = k2-1; 
                while (a(m(k2-1)) < b2) { ++k2; }
                auto r = parm/((1.0-parm)*fn); 
                auto is = 0; 
                auto sm = n-2*(k1-1);
                auto s = 0.0;
                for (int k = k1; k < k2; ++k) {
                    sm -= 2;
                    s = r * sm + am;
                    if (s > a(m(k-1)) && s <= a(m(k))) {
                        is = k;
                        break;
                    }
                }
                if (is) {
                    out = s;
                } else {
                    auto r2 = 2.0 * r; 
                    auto s1 = a(m(k1-1)); 
                    auto am2 = 2.0 * am;
                    auto cri = r2 * (a.array()-s1).abs().sum() + s1*(s1-am2); 
                    out = s1;
                    for (int k = k1+1; k < k2+1; ++k) {
                        s = a(m(k-1));
                        if (s == s1) continue;
                        auto c = r2 * (a.array()-s).abs().sum() + s*(s-am2);
                        if (c < cri) { cri = c; out = s; }
                        s1 = s;
                    }
                }
            }
        } 
    }
    out = std::max(
            (a.array()-cl(1)).maxCoeff(),
            std::min(
                (a.array()-cl(0)).minCoeff(), out) );
    return out;
}

template <class AType
        , class MType
        , class IntType
        , class ISType>
inline auto nintot(
        const AType& a,
        const MType& m,
        IntType nin,
        ISType&& is
        )
{
    auto nc = a.cols();
    is.setZero(); 
    int out = 0;
    for (int ic = 0; ic < nc; ++ic) {
        for (int j = 0; j < nin; ++j) {
            auto k = m(j)-1;
            if (is(k)) continue;
            if (!a(j,ic)) continue;
            is(k) = k+1;
            ++out;
        }
    }
    return out;
}

template <class FloatType
        , class ValueType
        , class XType
        , class YType
        , class GType
        , class WType
        , class JUType
        , class VPType
        , class CLType
        , class IntType
        , class ULamType
        , class A0Type
        , class AType
        , class MType
        , class KinType
        , class DevType
        , class ALMType>
inline void lognetn(
    ValueType beta,
    const XType& x,
    const YType& y,     // matrix 
    GType& g,           // matrix
    const WType& w,
    const JUType& ju,
    const VPType& vp,
    const CLType& cl,
    IntType ne,
    IntType nx,
    IntType nlam,
    ValueType flmin,
    const ULamType& ulam,
    ValueType shri, // basically thr
    bool isd,
    bool intr,
    IntType maxit,
    IntType kopt,
    IntType& lmu,
    A0Type& a0,     // matrix
    AType& a,       // matrix(nx, nc, nlam)
    MType& m,   // basically ia
    KinType& kin,
    ValueType& dev0,
    DevType& dev,
    ALMType& alm,
    IntType& nlp,
    IntType& jerr
        )
{
    // supplementary
    using int_t = IntType;
    using int_param_t = InternalParams;

    auto no = x.rows();
    auto ni = x.cols();
    auto nc = y.cols();
    // end supplementary

    // allocate inside elnet_point
    Eigen::MatrixXd b(ni+1, nc); b.setZero();
    Eigen::MatrixXd xv(ni, nc);
    Eigen::VectorXd ga(ni); ga.setZero();
    Eigen::MatrixXd bs(ni+1, nc); bs.setZero();
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXi is(std::max(nc, ni));
    Eigen::VectorXd sxp(no); sxp.setZero();
    Eigen::VectorXd sxpl(no);
    Eigen::VectorXd di(no);
    Eigen::VectorXi ixx(ni); ixx.setZero();
    Eigen::VectorXd r(no);
    Eigen::VectorXd v(no);
    Eigen::MatrixXd q(no, nc);
    // end allocate

    // initialize path config
    auto exmx = int_param_t::exmx;
    auto exmn = -int_param_t::exmx;
    auto pmin = int_param_t::pmin;
    auto pmax = 1.0 - pmin;
    auto emin = pmin / pmax;
    auto emax = 1.0/emin;
    auto pfm = (1.0 + pmin) * pmin;
    auto pfx = (1.0 - pmin) * pmax;
    auto vmin = pfm * pmax;
    auto bta = beta;
    auto omb = static_cast<FloatType>(1.0) - bta;
    auto dev1 = 0.0;
    dev0 = 0.0;

    // constructor of elnet_point
    for (int_t ic = 0; ic < nc; ++ic) {
        auto q0 = w.dot(y.col(ic));
        if (q0 <= pmin) {
            jerr = util::prob_min_reached_error(ic).err_code();
            return;
        }
        if (q0 >= 1.0 - pmin) {
            jerr = util::prob_max_reached_error(ic).err_code();
            return;
        }
        if (!intr) {
            q0 = 1.0 / nc;
            b(0, ic) = 0.0;
        }
        else {
            b(0, ic) = std::log(q0);
            dev1 -= q0 * b(0, ic);
        }
    }

    if (!intr) dev1 = std::log(nc);
    auto al = 0.0; // basically alm

    if ((g.array() == 0).all()) {
        b.row(0).array() -= b.row(0).sum() / nc;
        for (int_t ic = 0; ic < nc; ++ic) {
            q.col(ic).array() = std::exp(b(0, ic));
            sxp += q.col(ic);
        }
    } 
    else {
        for (int_t i = 0; i < no; ++i) {
            g.row(i).array() -= g.row(i).sum() / nc;
        }
        if (!intr) { 
            b.row(0).array() = 0.0; 
        } else {
            kazero(y, g, w, b.row(0));
        }
        dev1 = 0.0;
        for (int_t ic = 0; ic < nc; ++ic) {
            q.col(ic).array() = b(0,ic) + g.col(ic).array();
            dev1 -= w.dot( (y.col(ic).array() * q.col(ic).array()).matrix() );
            q.col(ic).array() = q.col(ic).array().exp();
            sxp += q.col(ic);
        }
        sxpl.array() = w.array() * sxp.array().log();
        for (int_t ic = 0; ic < nc; ++ic) {
            dev1 += y.col(ic).dot(sxpl);
        }
    }

    for (int_t ic = 0; ic < nc; ++ic) {
        for (int_t i = 0; i < no; ++i) {
            if (y(i,ic) > 0) dev0 += w(i) * y(i,ic) * std::log(y(i,ic));
        }
    }
    dev0 += dev1;

    if (kopt > 0) {
        if (isd > 0 && intr) { xv.array() = 0.25; }
        else { 
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j]) {
                    xv.row(j).array() = 0.25 * w.dot(x.col(j).array().square().matrix());
                }
            }
        }
    }

    // DEFINITELY initialize path config
    // "Begin: added by Naras"
    auto alf = 1.0;
    // "End: added by Naras"
    if (flmin < 1.0) { 
        auto eqs = std::max(int_param_t::eps,flmin); 
        alf = std::pow(eqs, 1.0/(nlam-1)); 
    }
    nlp = 0; 
    auto mnl = std::min(int_param_t::mnlam, nlam); 
    // end DEFINITELY
    
    auto shr = shri*dev0;
    m.setZero();    
    auto nin = 0;
    for (int_t ic = 0; ic < nc; ++ic) { 
        r.array() = w.array() * (y.col(ic).array() - q.col(ic).array() / sxp.array());
        for (int_t j = 0; j < ni; ++j) {
            if (!ju[j]) continue; 
            ga(j) = std::max(ga(j), std::abs(r.dot(x.col(j))));
        }
    }

    for (int_t ilm = 0; ilm < nlam; ++ilm) {
        // point config
        if (int_param_t::itrace) mock_setpb(ilm); 
        auto al0 = al;
        if (flmin >= 1.0) { al=ulam(ilm); }
        else if (ilm > 1) { al*=alf; }
        else if (ilm == 0) { al=int_param_t::big; }
        else { 
            al0 = 0.0;
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j] == 0) continue;
                if (vp(j) > 0.0) {
                    al0 = std::max(al0, ga(j) / vp(j));
                }
            }
            al0 /= std::max(bta,1.0e-3); 
            al = alf*al0;
        }
        auto al2 = al*omb; 
        auto al1 = al*bta; 
        // end point config
        
        Eigen::Map<Eigen::MatrixXd> a_slice(a.data() + nx * nc * ilm, nx, nc);
        auto tlam = bta*(2.0*al-al0);
        for (int_t k = 0; k < ni; ++k) {
            if (ixx(k) || !ju[k]) continue;
            if (ga(k) > tlam*vp(k)) ixx(k)=1;
        }
        bool ig = false;

        try {

            // :again:continue;
            while (1) {

                bool ix = false;
                bool jx = false;

                if (nlp > maxit) 
                    throw util::maxit_reached_error();
                
                try {
                    // iterate through each category
                    for (int_t ic = 0; ic < nc; ++ic) {

                        auto bs_ic = bs.col(ic);
                        auto b_ic = b.col(ic);
                        auto q_ic = q.col(ic);
                        auto y_ic = y.col(ic);
                        auto xv_ic = xv.col(ic);
                        auto g_ic = g.col(ic);

                        // SETUP FOR FITTING AT CATEGORY IC
                        bs_ic(0) = b_ic(0); 
                        if (nin > 0) {
                            for (int_t j = 0; j < nin; ++j) {
                                bs_ic(m(j)) = b_ic(m(j)); // should be correct since m contains 1-indexed values
                            }
                        }
                        auto xmz = 0.0;
                        for (int_t i = 0; i < no; ++i) {
                            auto pic = q_ic(i) / sxp(i);
                            if (pic < pfm) { pic = 0.0; v(i) = 0.0; }
                            else if (pic > pfx) { pic = 1.0; v(i) = 0.0; }
                            else { v(i) = w(i) * pic * (1.0 - pic); xmz += v(i); }
                            r(i) = w(i) * (y_ic(i) - pic);
                        }
                        if (xmz <= vmin) continue;
                        ig = true;
                        if (!kopt) {
                            for (int_t j = 0; j < ni; ++j) {
                                if (ixx(j)) xv_ic(j) = v.dot(x.col(j).array().square().matrix());
                            }
                        }
                        // END SETUP

                        // CD LOOP AT CATEGORY IC
                        while (1) {

                            ++nlp; 
                            auto dlx=0.0;
                            for (int_t k = 0; k < ni; ++k) {
                                if (!ixx(k)) continue;
                                auto bk = b_ic(k+1); 
                                auto gk = r.dot(x.col(k));
                                auto u = gk + xv_ic(k) * b_ic(k+1); 
                                auto au = std::abs(u)-vp(k)*al1;
                                if (au <= 0.0) { b_ic(k+1) = 0.0; }
                                else {
                                    b_ic(k+1) = std::max(cl(0,k),
                                            std::min(cl(1,k),
                                                std::copysign(au,u)/(xv_ic(k)+vp(k)*al2)));
                                }
                                // TODO: note that I changed the criterion here
                                // Potential issue with floating point precision errors.
                                if (bk == b_ic(k+1)) continue; 
                                auto d = b_ic(k+1)-bk; 
                                dlx = std::max(dlx, xv_ic(k)*d*d);
                                r.array() -= d * v.array() * x.col(k).array();
                                if (!mm(k)) { 
                                    ++nin; 
                                    if (nin > nx) throw util::max_active_reached_error();
                                    mm(k) = nin; 
                                    m(nin-1) = k+1;
                                }
                            }

                            auto d = 0.0; 
                            if (intr) d = r.sum()/xmz;
                            if (d) { 
                                b_ic(0) += d; 
                                dlx = std::max(dlx, xmz*d*d); 
                                r -= d*v;
                            }
                            if (dlx < shr) break; 
                            if (nlp > maxit) throw util::maxit_reached_error();
                            while (1) { 
                                ++nlp; 
                                auto dlx = 0.0;
                                for (int_t l = 0; l < nin; ++l) {
                                    auto k = m(l)-1;
                                    auto bk = b_ic(k+1);
                                    auto gk = r.dot(x.col(k));
                                    auto u = gk+xv_ic(k)*bk; 
                                    auto au = std::abs(u)-vp(k)*al1;
                                    if (au <= 0.0) { b_ic(k+1) = 0.0; }
                                    else {
                                      b_ic(k+1) = std::max(cl(0,k),
                                              std::min(cl(1,k),
                                                  std::copysign(au,u)/(xv_ic(k)+vp(k)*al2)));
                                    }
                                    if (b_ic(k+1) == bk) continue; 
                                    auto d = b_ic(k+1)-bk; 
                                    dlx = std::max(dlx,xv_ic(k)*d*d);
                                    r.array() -= d * v.array() * x.col(k).array();
                                }
                                auto d = 0.0; 
                                if (intr) d = r.sum()/xmz;
                                if (d) { 
                                    b_ic(0) += d; 
                                    dlx = std::max(dlx,xmz*d*d); 
                                    r -= d*v;
                                }
                                if (dlx < shr) break; 
                                if (nlp > maxit) throw util::maxit_reached_error();
                            } // end active cd-loop
                        } // end fitting for category ic

                        // TIDY-UP AFTER FITTING ALL CATEGORIES
                        auto d = b_ic(0) - bs_ic(0);
                        if (xmz * d * d > shr) ix = true;
                        if (!ix) {
                            for (int_t j = 0; j < nin; ++j) {
                                auto k = m(j)-1;
                                auto d = b_ic(k+1) - bs_ic(k+1);
                                if (xv_ic(k) * d * d > shr) { ix = true; break; }
                            }
                        }
                        for (int_t i = 0; i < no; ++i) {
                            auto fi = b_ic(0) + g_ic(i);
                            if (nin > 0) {
                                for (int_t j = 0; j < nin; ++j) {
                                    fi += b_ic(m(j)) * x(i, m(j)-1);
                                }
                            }
                            fi = std::min(std::max(exmn, fi), exmx);
                            sxp(i) -= q_ic(i);
                            q_ic(i) = std::min(std::max(emin * sxp(i), std::exp(fi)), emax*sxp(i));
                            sxp(i) += q_ic(i);
                        }
                        // END TIDY-UP
                    } // end loop over categories
                } // end try
                catch (const util::max_active_reached_error& e) {
                    jx = true;
                    --nin;  // necessary since otherwise the next part does 1 too many iterations.
                }

                auto s = -b.row(0).sum() / nc;
                b.row(0).array() += s;
                di.array() = s;
                for (int_t j = 0; j < nin; ++j) {
                    auto l = m(j)-1;
                    if (vp(l) <= 0) { s = b.row(l+1).sum()/nc; }
                    else { s = elc(beta, cl.col(l), b.row(l+1), is); }
                    b.row(l+1).array() -= s;
                    di -= s * x.col(l);
                }
                di.array() = di.array().exp();
                sxp.array() *= di.array();
                for (int_t ic = 0; ic < nc; ++ic) {
                    q.col(ic).array() *= di.array();
                }
                if (jx) throw util::max_active_reached_error();
                if (!ig) break;
                if (!ix) {
                    for (int_t k = 0; k < ni; ++k) {
                        if (ixx(k) || !ju[k]) continue;
                        ga(k) = 0.0;
                    }
                    for (int_t ic = 0; ic < nc; ++ic) {
                        r.array() = w.array() * (y.col(ic).array() - q.col(ic).array()/sxp.array());
                        for (int_t k = 0; k < ni; ++k) {
                            if (ixx(k) || !ju[k]) continue; 
                            ga(k) = std::max(ga(k), std::abs(r.dot(x.col(k))));
                        }
                    }
                    for (int_t k = 0; k < ni; ++k) {
                        if (ixx(k) || !ju[k]) continue;
                        if (ga(k) > al1*vp(k)) { 
                            ixx(k) = 1; 
                            ix = true;
                        }
                    }
                    if (!ix) break;
                }

            } // end outer-while :again:

        } // end try
        catch (const util::maxit_reached_error& e) {
            jerr = e.err_code(ilm);
            return;
        }
        catch (const util::elnet_error& e) {
            jerr = e.err_code(ilm);
            break;
        }
        
        // clean-up before next lambda iteration
        auto devi = 0.0;
        for (int_t ic = 0; ic < nc; ++ic) {
            if (nin > 0) {
                for (int_t j = 0; j < nin; ++j) {
                    a_slice(j,ic) = b(m(j), ic); 
                }
            }
            a0(ic, ilm) = b(0, ic);
            for (int_t i = 0; i < no; ++i) {
                if (y(i,ic) <= 0) continue;
                devi -= w(i) * y(i,ic) * std::log(q(i,ic) / sxp(i));
            }
        }
        kin(ilm) = nin;
        alm(ilm) = al; 
        lmu = ilm+1;
        dev(ilm) = (dev1-devi)/dev0; 
        if (!ig) break;
        if (lmu < mnl || flmin >= 1.0) continue; 
        auto prev_dev = (ilm == 0) ? 0 : dev(ilm-1);    // TODO: bug fix!
        if (nintot(a_slice,m,nin,is) > ne ||
            dev(ilm) > int_param_t::rsqmax ||
            dev(ilm) - prev_dev < int_param_t::sml) break;
        // end clean-up before next lambda iteration

    } // end for-loop path iteration

    // clean-up before finishing
    g.array() = q.array().log();
    for (int_t i = 0; i < no; ++i) {
        g.row(i).array() -= g.row(i).sum()/nc;
    }
    // end clean-up before finishing
}

} // namespace transl
} // namespace glmnetpp
