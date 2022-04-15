#pragma once
#include <Eigen/Core>
#include <testutil/translation/splognetn.hpp>
#include <testutil/translation/multelnet2.hpp> // chkbnds
#include <testutil/mock_pb.hpp>
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
        , class XVType
        , class XBType
        , class XSType
        , class A0Type
        , class AType
        , class MType
        , class KinType
        , class DevType
        , class ALMType>
inline void multsplognetn(
    ValueType beta,
    const XType& x,
    const YType& y,
    GType& g,
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
    bool intr,
    IntType maxit,
    const XVType& xv,
    const XBType& xb,
    const XSType& xs,
    IntType& lmu,
    A0Type& a0,
    AType& a,
    MType& m,   // basically ia
    KinType& kin,
    ValueType& dev0,
    DevType& dev,
    ALMType& alm,
    IntType& nlp,
    IntType& jerr
        )
{
    using int_param_t = InternalParams;

    auto ni = x.cols();
    auto nc = y.cols();
    auto no = x.rows();

    Eigen::MatrixXd b(ni+1, nc); b.setZero();
    Eigen::MatrixXd bs(ni+1, nc); bs.setZero();
    Eigen::MatrixXd q(no, nc);
    Eigen::MatrixXd r(no, nc);
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXd ga(ni); ga.setZero();
    Eigen::VectorXd gk(nc);
    Eigen::VectorXd del(nc);
    Eigen::VectorXd ixx(ni); ixx.setZero();
    Eigen::VectorXi is(std::max(nc, ni));
    Eigen::VectorXd sxp(no); sxp.setZero();
    Eigen::VectorXd sxpl(no);
    Eigen::VectorXd svr(nc); // TODO
    Eigen::VectorXd sc(no);  // TODO
    Eigen::VectorXd isc(nc);

    m.setZero();

    auto exmx = int_param_t::exmx;
    auto pmin = int_param_t::pmin;
    auto eps = int_param_t::eps;
    auto mnlam = int_param_t::mnlam;
    auto sml = int_param_t::sml;
    auto devmax = int_param_t::rsqmax;

    auto exmn = -exmx;
    auto pmax = 1.0-pmin; 
    auto emin = pmin/pmax; 
    auto emax = 1.0/emin;
    auto bta = beta; 
    auto omb=1.0-bta; 
    auto dev1=0.0; 
    dev0=0.0;
    for (int ic = 0; ic < nc; ++ic) {
        auto q0 = w.dot(y.col(ic));
        if (q0 <= pmin) { 
            jerr = util::prob_min_reached_error(ic).err_code(); 
            return; 
        }
        if (q0 >= pmax) {
            jerr = util::prob_max_reached_error(ic).err_code(); 
            return;
        } 
        if (!intr) { q0=1.0/nc; b(0,ic)=0.0; }
        else { b(0,ic)=std::log(q0); dev1-=q0*b(0,ic); }
    }

    if (!intr) dev1=std::log(nc); 
    if ((g.array() == 0.0).all()) {
       b.row(0).array() -= b.row(0).sum()/nc; 
       for (int ic = 0; ic < nc; ++ic) {
            q.col(ic).array() = std::exp(b(0,ic)); 
            sxp += q.col(ic);
       }
    }
    else {
        for (int i = 0; i < no; ++i) {
            g.row(i).array() -= g.row(i).sum() / nc;
        }
        if (intr) { 
            kazero(y,g,w,b.row(0)); 
        }
        dev1 = 0.0;
        for (int ic = 0; ic < nc; ++ic) {
            q.col(ic).array() = b(0,ic)+g.col(ic).array();
            dev1 -= w.dot((y.col(ic).array()*q.col(ic).array()).matrix());
            q.col(ic).array() = q.col(ic).array().exp(); 
            sxp += q.col(ic);
        }
        sxpl = w.array()*sxp.array().log(); 
        for (int ic = 0; ic < nc; ++ic) {
            dev1 += y.col(ic).dot(sxpl);
        }
    }
    for (int ic = 0; ic < nc; ++ic) {
        for (int i = 0; i < no; ++i) {
            if (y(i,ic) > 0.0) {
                dev0 += w(i)*y(i,ic)*std::log(y(i,ic));
            }
        }
    }
    dev0 += dev1;

    auto al=0.0;
    auto alf=1.0;
    if (flmin < 1.0) { 
        auto eqs=std::max(eps,flmin); 
        alf=std::pow(eqs, 1.0/(nlam-1)); 
    }
    auto nin=0; 
    nlp=0; 
    auto mnl=std::min(mnlam,nlam); 
    auto shr=shri*dev0;

    for (int ic = 0; ic < nc; ++ic) {
        r.col(ic) = w.array()*(y.col(ic).array()-q.col(ic).array()/sxp.array());
        svr(ic) = r.col(ic).sum(); // TODO
        for (int j = 0; j < ni; ++j) {
            if (ju[j]) {
                // TODO
                auto t = x.col(j).dot(r.col(ic));
                auto gk = (t-svr(ic)*xb(j))/xs(j);
                ga(j) += gk * gk;
            }
        }
    }
    if (apply_suggested_fixes) {
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            ga(j) = std::sqrt(ga(j));
        }
    } else {
        ga.array() = ga.array().sqrt();
    }

    for (int ilm = 0; ilm < nlam; ++ilm) {
        // point config
        if (int_param_t::itrace) mock_setpb(ilm); 
        auto al0 = al;
        if (flmin >= 1.0) { al=ulam(ilm); }
        else if (ilm > 1) { al*=alf; }
        else if (ilm == 0) { al=int_param_t::big; }
        else { 
            al0 = 0.0;
            for (int j = 0; j < ni; ++j) {
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
        for (int k = 0; k < ni; ++k) {
            if (ixx(k) || !ju[k]) continue;
            if (ga(k) > tlam*vp(k)) ixx(k)=1;
        }

        try {
            while (1) {
                bool ix = false;
                auto t = 0.0;
                if (nlp > maxit) { throw util::maxit_reached_error(); }

                for (int ic = 0; ic < nc; ++ic) {
                    t = std::max(t,
                            (q.col(ic).array()*
                             (1.0-q.col(ic).array()/sxp.array())/
                             sxp.array()).maxCoeff() );
                }
                if (t < eps) { 
                    throw util::below_min_variance_error();
                }
                t *= 2.0; 
                auto alt=al1/t; 
                auto al2t=al2/t;

                for (int ic = 0; ic < nc; ++ic) {
                    bs(0,ic)=b(0,ic); 
                    for (int j = 0; j < nin; ++j) {
                        bs(m(j),ic)=b(m(j),ic);
                    }

                    r.col(ic).array()=w.array()*(y.col(ic).array()-q.col(ic).array()/sxp.array())/t;
                    svr(ic) = r.col(ic).sum();

                    auto d=0.0; 
                    if(intr) d=svr(ic); // TODO
                    if (d) {
                       b(0,ic)+=d; 
                       r.col(ic)-=d*w; 
                       if (apply_suggested_fixes) {
                           svr(ic) = 0.0;
                       }
                    }
                }

                // BEGIN CD
                auto dlx = 0.0;
                while (1) {
                    nlp += nc; 
                    dlx=0.0;
                    for (int k = 0; k < ni; ++k) {
                        if (!ixx(k)) continue;
                        auto bk = b.row(k+1);
                        del.noalias() = bk.transpose();
                        // TODO:
                        gk = r.transpose() * x.col(k);
                        gk = (gk - svr * xb(k)) / xs(k);
                        gk += xv(k) * del;
                        auto gkn = gk.norm();
                        auto u = 1.0 - alt*vp(k)/gkn; 
                        if (u <= 0.0)  bk.setZero(); 
                        else {
                            bk = gk * (u/(xv(k)+vp(k)*al2t));
                            chkbnds(gk, gkn, xv(k), [&](auto i, auto) { return cl(i,k); },
                                     vp(k)*al2t, vp(k)*alt, bk, isc);
                        }
                        del = bk.transpose() - del;
                        if (del.array().abs().maxCoeff() <= 0.0) continue;
                        dlx = std::max(dlx, xv(k) * del.array().square().maxCoeff());
                        for (int ic = 0; ic < nc; ++ic) {
                            // TODO: current version is what I think should be the correct implementation.
                            // For documentation, we provide the direct translation of Fortran code.
                            // Testing shows this translation is correct.
                            if (apply_suggested_fixes) {
                                auto d_scaled = del(ic) / xs(k);
                                r.col(ic) -= d_scaled * x.col(k).cwiseProduct(w);
                                r.col(ic) += (d_scaled * xb(k)) * w;
                            } else {
                                auto jb = x.outerIndexPtr()[k];
                                auto je = x.outerIndexPtr()[k+1];
                                auto jx = x.innerIndexPtr();
                                for (; jb != je; ++jb) {
                                    r(jx[jb], ic) -= del(ic) * w(jx[jb]) * (x.valuePtr()[jb]-xb(k))/xs(k);
                                }
                            }
                        }
                        if (!mm(k)) { 
                            ++nin; 
                            if (nin > nx) throw util::max_active_reached_error();
                            mm(k) = nin; 
                            m(nin-1) = k+1;
                        }
                    }
                    
                    if (dlx < shr) break;
                    if (nlp > maxit) throw util::maxit_reached_error();

                    while (1) {
                        nlp += nc;
                        dlx = 0.0;
                        for (int l = 0; l < nin; ++l) {
                            auto k = m(l)-1;
                            auto bk = b.row(k+1);
                            del.noalias() = bk.transpose();

                            // TODO: optimize! product can go first
                            gk = r.transpose() * x.col(k);
                            gk = (gk - svr * xb(k)) / xs(k);
                            gk += xv(k) * del;
                            auto gkn = gk.norm();
                            auto u = 1.0 - alt*vp(k)/gkn; 
                            if (u <= 0.0)  bk.setZero(); 
                            else {
                                bk = gk * (u/(xv(k)+vp(k)*al2t));
                                chkbnds(gk, gkn, xv(k), [&](auto i, auto) { return cl(i,k); },
                                         vp(k)*al2t, vp(k)*alt, bk, isc);
                            }
                            del = bk.transpose() - del;
                            if (del.array().abs().maxCoeff() <= 0.0) continue;
                            dlx = std::max(dlx, xv(k) * del.array().square().maxCoeff());
                            for (int ic = 0; ic < nc; ++ic) {
                                // TODO: current version is what I think should be the correct implementation.
                                // For documentation, we provide the direct translation of Fortran code.
                                // Testing shows this translation is correct.
                                if (apply_suggested_fixes) {
                                    auto d_scaled = del(ic) / xs(k);
                                    r.col(ic) -= d_scaled * x.col(k).cwiseProduct(w);
                                    r.col(ic) += (d_scaled * xb(k)) * w;
                                } else {
                                    auto jb = x.outerIndexPtr()[k];
                                    auto je = x.outerIndexPtr()[k+1];
                                    auto jx = x.innerIndexPtr();
                                    for (; jb != je; ++jb) {
                                        r(jx[jb], ic) -= del(ic) * w(jx[jb]) * (x.valuePtr()[jb]-xb(k))/xs(k);
                                    }
                                }
                            }
                        }
                        if (dlx < shr) break;
                        if (nlp > maxit) throw util::maxit_reached_error();
                    }
                } // end cd

                for (int ic = 0; ic < nc; ++ic) {
                    auto d = b(0, ic) - bs(0, ic);
                    if (d * d > shr) ix = true;
                    if (!ix) {
                        for (int j = 0; j < nin; ++j) {
                            auto k = m(j)-1;
                            auto d = b(k+1,ic) - bs(k+1,ic);
                            if (xv(k) * d * d > shr) { ix = 1; break; }
                        }
                    }

                    // TODO
                    sc.array() = b(0,ic) + g.col(ic).array();
                    auto b0 = 0.0;
                    for (int j = 0; j < nin; ++j) {
                        auto l = m(j)-1;
                        auto d_scaled = b(l+1,ic) / xs(l);
                        sc += d_scaled * x.col(l);
                        b0 -= d_scaled * xb(l);
                    }
                    sc.array() = (sc.array()+b0).max(exmn).min(exmx);
                    sxp -= q.col(ic);
                    q.col(ic).array() = (emin*sxp.array()).max(sc.array().exp()).min(emax*sxp.array());
                    sxp += q.col(ic);
                }
                auto s = -b.row(0).sum()/nc;
                b.row(0).array() += s;
                if (!ix) {
                    for (int k = 0; k < ni; ++k) {
                        if (ixx(k) || !ju[k]) continue;
                        ga(k) = 0.0;
                    } 
                    for (int ic = 0; ic < nc; ++ic) {
                        r.col(ic) = w.cwiseProduct(y.col(ic) - q.col(ic).cwiseQuotient(sxp));
                        if (apply_suggested_fixes) {
                            svr(ic) = r.col(ic).sum();
                        }
                        for (int k = 0; k < ni; ++k) {
                            if (ixx(k) || !ju[k]) continue;
                            // TODO
                            auto gk = x.col(k).dot(r.col(ic));
                            gk = (gk - svr(ic) * xb(k))/xs(k);
                            ga(k) += gk * gk;
                        }
                    }
                    if (apply_suggested_fixes) {
                        for (int j = 0; j < ni; ++j) {
                            if (ixx(j) || !ju[j]) continue;
                            ga(j) = std::sqrt(ga(j));
                        }
                    } else {
                        ga.array() = ga.array().sqrt();
                    }
                    for (int k = 0; k < ni; ++k) {
                        if (ixx(k) || !ju[k]) continue;
                        if (ga(k) > al1 * vp(k)) {
                            ixx(k) = 1;
                            ix = true;
                        }
                    }
                    if (ix) continue;
                    break;
                }

            } // end irls loop

        } // end try
        catch (const util::maxit_reached_error& e) {
            jerr = e.err_code(ilm);
            return;
        } 
        catch (const util::bnorm_maxit_reached_error& e) {
            jerr = e.err_code(ilm);
            return;
        }
        catch (const util::max_active_reached_error& e) {
            jerr = e.err_code(ilm);
            break;
        } 
        catch (const util::below_min_variance_error& e) {
            jerr = e.err_code(ilm);
            break;
        }

        auto devi = 0.0;
        for (int ic = 0; ic < nc; ++ic) {
            for (int j = 0; j < nin; ++j) {
                a_slice(j, ic) = b(m(j), ic);
            }
            a0(ic, ilm) = b(0, ic);
            for (int i = 0; i < no; ++i) {
                if (y(i,ic) <= 0.0) continue;
                devi -= w(i) * y(i,ic) * std::log(q(i,ic)/sxp(i));
            }
        }
        kin(ilm)=nin;
        alm(ilm)=al;
        lmu=ilm+1;
        dev(ilm) = (dev1-devi)/dev0;
        if (lmu < mnl || flmin >= 1.0) continue;
        auto me = (a_slice.col(0).head(nin).array() != 0.0).count();
        auto prev_dev = (ilm <= 0) ? 0.0 : dev(ilm-1);
        if (me > ne ||
            dev(ilm) > devmax ||
            dev(ilm) - prev_dev < sml) break;

    } // end for-loop over lambdas

    g = q.array().log().matrix();
    for (int i = 0; i < no; ++i) {
        g.row(i).array() -= g.row(i).sum()/nc;
    }
}
    
} // namespace transl
} // namespace glmnetpp
