#pragma once
#include <Eigen/Core>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/lognetn.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>

namespace glmnetpp {
namespace transl {

template <class FloatType
        , bool apply_suggested_fixes
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
        , class XBType
        , class XSType
        , class A0Type
        , class AType
        , class MType
        , class KinType
        , class DevType
        , class ALMType>
inline void splognetn(
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
    const XBType& xb, // TODO: mean of columns of x
    const XSType& xs, // TODO: std of columns of x
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
    Eigen::MatrixXd bs(ni+1, nc); bs.setZero();
    Eigen::MatrixXd q(no, nc);
    Eigen::VectorXd xm(ni+1); // TODO
    Eigen::VectorXd r(no);
    Eigen::VectorXd v(no);
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXd ga(ni); ga.setZero();
    Eigen::VectorXi ixx(ni); ixx.setZero();
    Eigen::VectorXi is(std::max(nc, ni));
    Eigen::VectorXd sxp(no); sxp.setZero();
    Eigen::VectorXd sxpl(no);
    Eigen::VectorXd di(no);
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
                    // TODO
                    auto xj_sq = x.col(j).cwiseProduct(x.col(j));
                    if (apply_suggested_fixes) {
                        xv.row(j).array() = 0.25 * 
                            (xj_sq.dot(w) - xb(j) * xb(j))/(xs(j)*xs(j));
                    } else {
                        xv.row(j).array() = 0.25 * 
                            (xj_sq.dot(w) - xb(j) * xb(j));
                    } 
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
    auto svr = 0.0; // TODO
    auto o = 0.0; // TODO
    for (int_t ic = 0; ic < nc; ++ic) { 
        v.array() = q.col(ic).array()/sxp.array(); // TODO
        r.array() = w.array() * (y.col(ic).array() - v.array()); // TODO
        v.array() = w.array() * v.array() * (1.0 - v.array()); // TODO
        if (apply_suggested_fixes) {
            svr = r.sum();
        }
        for (int_t j = 0; j < ni; ++j) {
            if (!ju[j]) continue; 
            // TODO
            auto grad = x.col(j).dot(
                    (r.array() + v.array() * o).matrix()
                    );
            ga(j) = std::max(ga(j), std::abs((grad-svr*xb(j))/xs(j)));
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
                        xm(0) = 0.0; // TODO
                        svr = 0.0; // TODO
                        o = 0.0; // TODO
                        for (int_t i = 0; i < no; ++i) {
                            auto pic = q_ic(i) / sxp(i);
                            if (pic < pfm) { pic = 0.0; v(i) = 0.0; }
                            else if (pic > pfx) { pic = 1.0; v(i) = 0.0; }
                            else { v(i) = w(i) * pic * (1.0 - pic); xm(0) += v(i); } // TODO
                            r(i) = w(i) * (y_ic(i) - pic);
                            svr += r(i); // TODO
                        }
                        if (xm(0) <= vmin) continue; // TODO
                        ig = true;

                        // TODO
                        for (int_t j = 0; j < ni; ++j) {
                            if (!ixx(j)) continue;
                            xm(j+1) = x.col(j).dot(v);
                            if (!kopt) {
                                auto xj_sq = x.col(j).cwiseProduct(x.col(j));
                                xv_ic(j) = (xj_sq.dot(v) - 2.0*xb(j)*xm(j+1)+xm(0)*xb(j)*xb(j))/(xs(j) * xs(j));
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
                                // TODO
                                auto gk = x.col(k).dot((r.array() + v.array() * o).matrix());
                                gk = (gk - svr * xb(k))/xs(k);

                                auto u = gk + xv_ic(k) * b_ic(k+1); 
                                auto au = std::abs(u)-vp(k)*al1;
                                if (au <= 0.0) { b_ic(k+1) = 0.0; }
                                else {
                                    b_ic(k+1) = std::max(cl(0,k),
                                            std::min(cl(1,k),
                                                std::copysign(au,u)/(xv_ic(k)+vp(k)*al2)));
                                }
                                if (bk == b_ic(k+1)) continue; 
                                auto d = b_ic(k+1)-bk; 
                                dlx = std::max(dlx, xv_ic(k)*d*d);
                                if (!mm(k)) { 
                                    ++nin; 
                                    if (nin > nx) throw util::max_active_reached_error();
                                    mm(k) = nin; 
                                    m(nin-1) = k+1;
                                    xm(k+1) = x.col(k).dot(v); // TODO
                                }
                                // TODO
                                auto d_scaled = d / xs(k);
                                r.array() -= d_scaled * x.col(k).cwiseProduct(v);
                                o += d_scaled * xb(k);
                                svr -= d_scaled * (xm(k+1) - xb(k) * xm(0));
                            }

                            auto d = 0.0; 
                            if (intr) d = svr/xm(0); // TODO
                            if (d) { 
                                b_ic(0) += d; 
                                dlx = std::max(dlx, xm(0)*d*d); // TODO
                                r -= d*v;
                                if (apply_suggested_fixes) {
                                    svr = 0.0; 
                                } else {
                                    svr -= d * xm(0); // TODO
                                }
                            }
                            if (dlx < shr) break; 
                            if (nlp > maxit) throw util::maxit_reached_error();
                            while (1) { 
                                ++nlp; 
                                auto dlx = 0.0;
                                for (int_t l = 0; l < nin; ++l) {
                                    auto k = m(l)-1;
                                    auto bk = b_ic(k+1);
                                    // TODO
                                    auto gk = x.col(k).dot((r.array() + v.array()*o).matrix());
                                    gk = (gk-svr*xb(k))/xs(k);

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

                                    // TODO
                                    auto d_scaled = d / xs(k);
                                    r.array() -= d_scaled * x.col(k).cwiseProduct(v);
                                    o += d_scaled * xb(k);
                                    svr -= d_scaled * (xm(k+1) - xb(k) * xm(0));
                                }
                                auto d = 0.0; 
                                if (intr) d = svr/xm(0); // TODO
                                if (d) { 
                                    b_ic(0) += d; 
                                    dlx = std::max(dlx,xm(0)*d*d); // TODO
                                    r -= d*v;
                                    if (apply_suggested_fixes) {
                                        svr = 0.0;
                                    } else {
                                        svr -= d * xm(0); // TODO
                                    }
                                }
                                if (dlx < shr) break; 
                                if (nlp > maxit) throw util::maxit_reached_error();
                            } // end active cd-loop
                        } // end fitting for category ic

                        // TIDY-UP AFTER FITTING ALL CATEGORIES
                        auto d = b_ic(0) - bs_ic(0);
                        if (xm(0) * d * d > shr) ix = true; // TODO
                        if (!ix) {
                            for (int_t j = 0; j < nin; ++j) {
                                auto k = m(j)-1;
                                auto d = b_ic(k+1) - bs_ic(k+1);
                                if (xv_ic(k) * d * d > shr) { ix = true; break; }
                            }
                        }

                        // TODO
                        di.array() = b_ic(0) + g_ic.array();
                        auto b0 = 0.0;
                        for (int_t j = 0; j < nin; ++j) {
                            auto l = m(j)-1;
                            auto b_scaled = b_ic(l+1) / xs(l);
                            di += b_scaled * x.col(l);
                            b0 -= b_scaled * xb(l);
                        }
                        di.array() = (di.array()+b0).max(exmn).min(exmx);
                        sxp -= q.col(ic);
                        q.col(ic).array() = (emin*sxp.array()).max(di.array().exp()).min(emax*sxp.array());
                        sxp += q.col(ic);
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
                auto b0 = 0.0; // TODO
                for (int_t j = 0; j < nin; ++j) {
                    auto l = m(j)-1;
                    if (vp(l) <= 0) { s = b.row(l+1).sum()/nc; }
                    else { s = elc(beta, cl.col(l), b.row(l+1), is); }
                    b.row(l+1).array() -= s;
                    // TODO
                    auto s_scaled = s/xs(l);
                    di -= s_scaled * x.col(l);
                    b0 += s_scaled * xb(l);
                }
                di.array() = (di.array() + b0).exp(); // TODO
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
                        // TODO
                        v.array() = q.col(ic).array() / sxp.array(); 
                        r.array() = w.array() * (y.col(ic).array() - v.array());
                        v.array() = w.array() * v.array() * (1.0-v.array());
                        if (apply_suggested_fixes) {
                            svr = r.sum();
                        }
                        for (int_t k = 0; k < ni; ++k) {
                            if (ixx(k) || !ju[k]) continue; 
                            // TODO
                            auto grad = x.col(k).dot(
                                    (r.array() + v.array()*o).matrix()
                                    );
                            ga(k) = std::max(ga(k), std::abs((grad-svr*xb(k))/xs(k)) );
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
