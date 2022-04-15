#pragma once
#include <Eigen/Core>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/lognet2n.hpp>
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
inline void splognet2n(
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
    bool isd,
    bool intr,
    IntType maxit,
    IntType kopt,
    const XBType& xb, // I think column means of X?
    const XSType& xs, // I think column stddevs of X?
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
    using int_t = IntType;
    using int_param_t = InternalParams;

    auto no = x.rows();
    auto ni = x.cols();

    Eigen::VectorXd b(ni+1); b.setZero();
    Eigen::VectorXd xm(ni+1);
    Eigen::VectorXd xv(ni);
    Eigen::VectorXd bs(ni+1); bs.setZero();
    Eigen::VectorXd ga(ni);
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXi ixx(ni); ixx.setZero();
    Eigen::VectorXd q(no);
    Eigen::VectorXd r(no);
    Eigen::VectorXd v(no);
    Eigen::VectorXd sc(no);

    auto fmax = std::log(
            static_cast<FloatType>(1.0)/int_param_t::pmin 
            - static_cast<FloatType>(1.0));
    auto fmin = -fmax;
    auto vmin = (static_cast<FloatType>(1.0) + int_param_t::pmin) 
        * int_param_t::pmin 
        * (static_cast<FloatType>(1.0) - int_param_t::pmin);
    auto bta = beta;
    auto omb = static_cast<FloatType>(1.0) - bta;
    auto q0 = w.dot(y);
    if (q0 <= int_param_t::pmin) {
        jerr = util::prob_min_reached_error(0).err_code();
        return;
    }
    if (q0 >= static_cast<FloatType>(1.0)-int_param_t::pmin) {
        jerr = util::prob_max_reached_error(0).err_code();
        return;
    }

    if (!intr) q0 = static_cast<FloatType>(0.5);
    auto bz = 0.0; 
    if (intr) bz = std::log(q0/(static_cast<FloatType>(1.0)-q0));

    auto dev1 = 0.0;
    if ((g.array() == 0).all()) {
        auto vi = q0*(static_cast<FloatType>(1.0)-q0); 
        b(0) = bz; 
        v = vi * w;
        r.array() = w.array() * (y.array()-q0); 
        q.array() = q0; 
        xm(0) = vi;  // TODO
        dev1 = -(bz*q0 + std::log(static_cast<FloatType>(1.0)-q0));
    } 

    else {
        b(0) = 0.0;
        if (intr) { 
            b(0) = azero(y,g,w); 
        }
        q.array() = 1.0 / (1.0 + (-b(0)-g.array()).exp()); 
        v.array() = w.array() * q.array() * (1.0-q.array()); 
        r.array() = w.array() * (y-q).array(); 
        xm(0) = v.sum(); // TODO
        dev1 = -(b(0)*q0 + 
                w.dot( (y.array()*g.array() + 
                        (static_cast<FloatType>(1.0)-q.array()).log()).matrix() ));
    }

    if (kopt > 0) {
        if (isd > 0 && intr) { xv.array() = 0.25; }
        else { 
            for (int_t j = 0; j < ni; ++j) {
                if (ju[j]) {
                    // TODO
                    auto xj_sq = x.col(j).cwiseProduct(x.col(j));
                    if (apply_suggested_fixes) {
                        xv(j) = 0.25 *
                            (xj_sq.dot(w) - xb(j) * xb(j))/(xs(j)*xs(j));
                    } else {
                        xv(j) = 0.25 * 
                            (xj_sq.dot(w) - xb(j) * xb(j));
                    }
                }
            }
        }
    }
    
    dev0 = dev1;
    for (int_t i = 0; i < no; ++i) {
        if (y(i) > 0.0) dev0 += w(i)*y(i)*std::log(y(i));
        else if (y(i) < 1.0) {
            dev0 += w(i)*(static_cast<FloatType>(1.0)-y(i))
                    *std::log(static_cast<FloatType>(1.0)-y(i));
        }
    }
    // "Begin: added by Naras"
    auto alf = 1.0;
    // "End: added by Naras"
    if (flmin < 1.0) { 
        auto eqs = std::max(int_param_t::eps,flmin); 
        alf = std::pow(eqs, 1.0/(nlam-1)); 
    }
    m.setZero();    
    auto nin = 0;
    auto o = 0.0; // TODO
    auto svr = 0.0; // TODO
    if (apply_suggested_fixes){
        svr = r.sum();
    }
    nlp = 0; 
    auto mnl = std::min(int_param_t::mnlam, nlam); 
    auto shr = shri*dev0;
    auto al = 0.0; // basically alm
    for (int_t j = 0; j < ni; ++j) {
        if (!ju[j]) continue; 
        // TODO
        auto grad = x.col(j).dot(
                (r.array() + v.array() * o).matrix()
                );
        ga(j) = std::abs((grad - svr*xb(j)) / xs(j));
    }

    for (int_t ilm = 0; ilm < nlam; ++ilm) {
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
        auto tlam = bta*(2.0*al-al0);
        for (int_t k = 0; k < ni; ++k) {
            if (ixx(k) || !ju[k]) continue;
            if (ga(k) > tlam*vp(k)) ixx(k)=1;
        }

        try {

            // :again:continue;
            while (1) {

                if (nlp > maxit) 
                    throw util::maxit_reached_error();
                
                // update bs with b
                bs(0) = b(0); 
                if (nin > 0) {
                    for (int_t j = 0; j < nin; ++j) {
                        bs(m(j)) = b(m(j)); // should be correct since m contains 1-indexed values
                    }
                }

                for (int_t j = 0; j < ni; ++j) {
                    if (!ixx(j)) continue;
                    // TODO: supposed to put into sc
                    xm(j+1)=x.col(j).dot(v);
                    if (!kopt) {
                        auto xj_sq = x.col(j).cwiseProduct(x.col(j));
                        xv(j) = xj_sq.dot(v);
                        xv(j)=(xv(j)-2.0*xb(j)*xm(j+1)+xm(0)*xb(j)*xb(j))/(xs(j) * xs(j));
                    }
                }

                while (1) {

                    ++nlp; 
                    auto dlx=0.0;
                    for (int_t k = 0; k < ni; ++k) {
                        if (!ixx(k)) continue;
                        auto bk = b(k+1); 
                        // TODO
                        auto gk = x.col(k).dot(
                                (r.array() + v.array() * o).matrix()
                                );
                        gk = (gk - svr * xb(k)) / xs(k);

                        auto u = gk + xv(k) * bk; 
                        auto au = std::abs(u)-vp(k)*al1;
                        if (au <= 0.0) { b(k+1) = 0.0; }
                        else {
                            b(k+1) = std::max(cl(0,k),
                                    std::min(cl(1,k),
                                        std::copysign(au,u)/(xv(k)+vp(k)*al2)));
                        }
                        if (bk == b(k+1)) continue; 
                        auto d = b(k+1)-bk; 
                        dlx = std::max(dlx, xv(k)*d*d);
                        if (!mm(k)) { 
                            ++nin; 
                            if (nin > nx) throw util::max_active_reached_error();
                            mm(k) = nin; 
                            m(nin-1) = k+1;
                            xm(k+1) = x.col(k).dot(v); // TODO
                        }
                        // TODO
                        auto d_scaled = d / xs(k);
                        r -= d_scaled * x.col(k).cwiseProduct(v);
                        o += d_scaled * xb(k); 
                        svr -= d_scaled*(xm(k+1)-xb(k)*xm(0));
                    }

                    auto d = 0.0; 
                    if (intr) d = svr/xm(0); // TODO
                    if (d) { 
                        b(0) += d; 
                        dlx = std::max(dlx, xm(0)*d*d);  // TODO
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
                            auto bk = b(k+1);
                            // TODO
                            auto gk = x.col(k).dot(
                                    (r.array() + v.array() * o).matrix()
                                    );
                            gk = (gk - svr * xb(k)) / xs(k);

                            auto u = gk+xv(k)*bk; 
                            auto au = std::abs(u)-vp(k)*al1;
                            if (au <= 0.0) { b(k+1) = 0.0; }
                            else {
                              b(k+1) = std::max(cl(0,k),
                                      std::min(cl(1,k),
                                          std::copysign(au,u)/(xv(k)+vp(k)*al2)));
                            }
                            auto d = b(k+1)-bk; 
                            if (b(k+1) == bk) continue; 
                            dlx = std::max(dlx,xv(k)*d*d);
                            // TODO
                            auto d_scaled = d / xs(k);
                            r -= d_scaled * x.col(k).cwiseProduct(v);
                            o += d_scaled * xb(k); 
                            svr -= d_scaled*(xm(k+1)-xb(k)*xm(0));
                        }
                        auto d = 0.0; 
                        if (intr) d = svr/xm(0); // TODO
                        if (d) { 
                            b(0) += d; 
                            dlx = std::max(dlx,xm(0)*d*d); // TODO
                            r -= d*v;
                            if (apply_suggested_fixes) {
                                svr = 0.0;
                            } else {
                                svr -= d * xm(0);
                            }
                        }
                        if (dlx < shr) break; 
                        if (nlp > maxit) throw util::maxit_reached_error();
                    }
                }

                sc.array() = b(0); 
                auto b0=0.0;
                for (int j = 0; j < nin; ++j) {
                    auto l=m(j)-1; 
                    sc += (b(l+1) / xs(l)) * x.col(l);
                    b0 -= b(l+1)*xb(l)/xs(l);
                }
                sc.array() += b0;

                for (int_t i = 0; i < no; ++i) {
                    auto fi = sc(i)+g(i); // TODO
                    if (fi < fmin) { q(i)=0.0; } 
                    else if (fi > fmax) { q(i) = 1.0; }
                    else { q(i)=1.0/(1.0+std::exp(-fi)); }
                }
                v.array() = w.array()*q.array()*(1.0-q.array()); 
                xm(0) = v.sum(); // TODO
                if (xm(0) < vmin) break; // TODO
                r.array() = w.array() * (y-q).array();
                svr = r.sum(); // TODO
                o = 0.0; // TODO
                auto diff0 = b(0) - bs(0);
                if (xm(0) * diff0 * diff0 < shr) { // TODO
                    bool ix = false;
                    for (int_t j = 0; j < nin; ++j) {
                        auto k = m(j)-1;
                        auto diff = b(k+1) - bs(k+1);
                        if (xv(k)*diff*diff < shr) continue; 
                        ix = true; 
                        break;
                    }
                    if (!ix) {
                        for (int_t k = 0; k < ni; ++k) {
                            if (ixx(k) || !ju[k]) continue; 
                            // TODO:
                            auto grad = x.col(k).dot(
                                    (r.array() + v.array() * o).matrix()
                                    );
                            ga(k) = std::abs((grad - svr * xb(k)) / xs(k));

                            if (ga(k) > al1*vp(k)) { 
                                ixx(k) = 1; 
                                ix = true;
                            }
                        }
                        if (!ix) break;
                    }
                } // end if-block if converged

            } // end outer-while :again

        } // end try
        catch (const util::maxit_reached_error& e) {
            jerr = e.err_code(ilm);
            return;
        }
        catch (const util::elnet_error& e) {
            jerr = e.err_code(ilm);
            break;
        }

        if (nin > 0) {
            for (int_t j = 0; j < nin; ++j) {
                a(j,ilm) = b(m(j)); 
            }
        }
        kin(ilm) = nin;
        a0(ilm) = b(0); 
        alm(ilm) = al; 
        lmu = ilm+1;
        auto devi = dev2<FloatType>(w,y,q,int_param_t::pmin);
        dev(ilm) = (dev1-devi)/dev0; 
        if (lmu < mnl || flmin >= 1.0) continue; 
        auto me = (a.col(ilm).head(nin).array() != 0.0).count();
        auto prev_dev = (ilm == 0) ? 0 : dev(ilm-1);    // TODO: bug fix!
        if (me > ne ||
            dev(ilm) > int_param_t::rsqmax ||
            dev(ilm) - prev_dev < int_param_t::sml) break;
        if (xm(0) < vmin) break; // TODO

    } // end for-loop path iteration

    g.array() = (q.array()/(1.0-q.array())).log();
}

} // namespace transl
} // namespace glmnetpp
