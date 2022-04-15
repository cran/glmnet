#pragma once
#include <limits>
#include <testutil/mock_pb.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <testutil/internal.hpp>
#include <Eigen/Core>

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
        , class CAType
        , class IAType
        , class KinType
        , class DevType
        , class AlmType>
inline void spfishnet1(
        ValueType beta,
        const XType& x,
        const YType& y,
        GType& g,
        const WType& q,
        const JUType& ju,
        const VPType& vp,
        const CLType& cl,
        IntType ne,
        IntType nx,
        IntType nlam,
        ValueType flmin,
        const ULamType& ulam,
        ValueType thr,
        bool intr,
        IntType maxit,
        const XBType& xb,
        const XSType& xs,
        IntType& lmu,
        A0Type& a0,
        CAType& ca,
        IAType& ia,
        KinType& kin,
        ValueType& dev0,
        DevType& dev,
        AlmType& alm,
        IntType& nlp,
        IntType& jerr
        )
{
    using int_param_t = InternalParams;    

    auto no = x.rows();
    auto ni = x.cols();
    auto sml = int_param_t::sml;
    auto eps = int_param_t::eps;
    auto mnlam = int_param_t::mnlam;
    auto itrace = int_param_t::itrace;
    auto devmax = int_param_t::rsqmax;
    sml *= 10.0;

    Eigen::VectorXd a(ni); a.setZero();
    Eigen::VectorXd as(ni); as.setZero();
    Eigen::VectorXd t(no);
    Eigen::VectorXi mm(ni); mm.setZero();
    Eigen::VectorXd ga(ni);
    Eigen::VectorXi ixx(ni); ixx.setZero();
    Eigen::VectorXd wr(no);
    Eigen::VectorXd v(ni);
    Eigen::VectorXd w(no);
    // Eigen::VectorXd f(no); // TODO: remove
    Eigen::VectorXd qy(no); // TODO
    Eigen::VectorXd xm(ni); // TODO

    auto bta = beta; 
    auto omb = 1.0 - bta;

    // TODO
    //t.array()=q.array()*y.array(); 
    //auto yb=t.sum(); 
    qy.array()=q.array()*y.array(); 
    auto yb=qy.sum(); 

    auto fmax=std::log(std::numeric_limits<double>::max()*0.1);
    double dv0 = 0;
    double az = 0;
    double uu = 0; // TODO
    if ((g.array() == 0).all()) {
        t.setZero();
        if (intr) { 
            w=q*yb;  
            az=std::log(yb); 
            uu = az; // TODO
            xm = yb * xb; // TODO
            // f.array()=az; 
            dv0=yb*(az-static_cast<FloatType>(1.0)); 
        }
        else { 
            w=q; 
            xm.setZero(); // TODO
            az = 0;
            uu = 0; // TODO
            //f.setZero(); 
            dv0=-static_cast<FloatType>(1.0);
        }
    }
    else {
        w.array() = q.array() * (g.array().abs().min(fmax)).matrix().binaryExpr(g,
                    [&](auto x, auto y) { return std::copysign(x, y); }).array().exp();
        auto v0=w.sum();
        t = g;
        if (intr) { 
            auto eaz=yb/v0; 
            w *= eaz; 
            az=std::log(eaz); 
            uu = az; // TODO
            //f.array()=az+g.array();
            //dv0=t.dot(g)-yb*(1.0-az);
            dv0=qy.dot(g)-yb*(1.0-az);
        }
        else { 
            az=0.0; 
            uu = 0.0; // TODO
            //f=g; 
            //dv0=t.dot(g)-v0; 
            dv0=qy.dot(g)-v0; 
        }
        // TODO
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = x.col(j).dot(w);
        }
    }

    auto v0=1.0; 
    if (intr) v0=yb; 
    double tt = 0.0;
    if (apply_suggested_fixes) {
        tt = yb - v0 * (static_cast<FloatType>(1.0)-uu);
    } else {
        tt = yb * uu; // TODO
    }
    //wr=t-w; 
    if (apply_suggested_fixes) {
        wr = qy - w * (static_cast<FloatType>(1.0)-uu);
    } else {
        wr = qy - q * (yb * (static_cast<FloatType>(1.0)-uu)); // TODO
    }
    auto dvr=-yb;
    for (int i = 0; i < no; ++i) {
        // TODO
        //if (t(i) > 0.0) dvr += t(i)*std::log(y(i));
        if (qy(i) > 0.0) dvr += qy(i)*std::log(y(i));
    }
    dvr-=dv0; 
    dev0=dvr;
    auto alf=1.0;
    if (flmin < 1.0) { 
        auto eqs=std::max(eps,flmin); 
        alf=std::pow(eqs, (1.0/(nlam-1))); 
    }
    ia.setZero(); 
    nlp = 0;
    auto nin = 0; 
    auto mnl=std::min(mnlam,nlam); 
    thr*=dev0; 
    auto al=0.0;
    for (int j = 0; j < ni; ++j) {
        if (!ju[j]) continue; 
        // TODO
        //ga(j)=std::abs(wr.dot(x.col(j)));
        ga(j) = std::abs(x.col(j).dot(wr) - uu*(xm(j)-v0*xb(j)) - xb(j)*tt) / xs(j);
    }

    for (int ilm = 0; ilm < nlam; ++ilm) {
        if (itrace) mock_setpb(ilm); 
        auto al0=al;
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
        auto tlam = bta*(static_cast<FloatType>(2.0)*al-al0);
        for (int k = 0; k < ni; ++k) {
            if (ixx(k) || !ju[k]) continue;
            if (ga(k) > tlam*vp(k)) ixx(k)=1;
        }

        try {
            while (1) {
                if (nlp > maxit) throw util::maxit_reached_error();
                auto az0=az; 
                for (int i = 0; i < nin; ++i) {
                    as(ia(i)-1) = a(ia(i)-1);
                }
                for (int j = 0; j < ni; ++j) {
                    if (ixx(j)) {
                        // TODO
                        // v(j) = w.dot(x.col(j).array().square().matrix());
                        xm(j) = x.col(j).dot(w);
                        v(j) = (x.col(j).cwiseProduct(x.col(j)).dot(w) 
                                - static_cast<FloatType>(2.0)*xb(j)*xm(j) 
                                + v0*xb(j)*xb(j))/(xs(j)*xs(j));
                    }
                }

                // CD loop
                while (1) {
                    ++nlp; 
                    auto dlx=0.0;
                    for (int k = 0; k < ni; ++k) {
                        if (!ixx(k)) continue; 
                        auto ak=a(k);
                        // TODO
                        auto gk = (x.col(k).dot(wr) - uu*(xm(k)-v0*xb(k))-xb(k)*tt)/xs(k);
                        auto u=gk+v(k)*ak; 
                        auto au=std::abs(u)-vp(k)*al1;
                        if (au <= 0.0) { a(k)=0.0; }
                        else {
                           a(k)=std::max(cl(0,k),
                                   std::min(cl(1,k),
                                       std::copysign(au,u)/(v(k)+vp(k)*al2)));
                        }
                        if (a(k) == ak) continue; 
                        if (!mm(k)) { 
                            ++nin; 
                            if (nin > nx) throw util::max_active_reached_error();
                            mm(k)=nin; 
                            ia(nin-1)=k+1;
                        }
                        auto d=a(k)-ak; 
                        dlx=std::max(dlx,v(k)*d*d);

                        // TODO
                        auto d_scaled = d/xs(k);
                        wr -= d_scaled*x.col(k).cwiseProduct(w);
                        uu -= d_scaled * xb(k);
                        tt -= d_scaled * xm(k);
                        //wr.array()-=d*w.array()*x.col(k).array(); 
                        //f += d*x.col(k);
                    }
                    if (intr) { 
                        // TODO
                        //auto d=wr.sum()/v0;
                        auto d = tt/v0 - uu;

                        az += d; 
                        dlx = std::max(dlx,v0*d*d); 
                        uu += d;
                        //wr-=d*w; 
                        //f.array()+=d;
                    }
                    if (dlx < thr) break; 
                    if (nlp > maxit) throw util::maxit_reached_error(); 
                    while (1) { 
                        ++nlp; 
                        auto dlx=0.0;
                        for (int l = 0; l < nin; ++l) {
                            auto k=ia(l)-1; 
                            auto ak=a(k);
                            // TODO
                            // auto gk = wr.dot(x.col(k));
                            auto gk = (x.col(k).dot(wr) - uu*(xm(k)-v0*xb(k))-xb(k)*tt)/xs(k);
                            auto u=gk+v(k)*ak; 
                            auto au=std::abs(u)-vp(k)*al1;
                            if (au <= 0.0) { a(k)=0.0; }
                            else {
                                a(k)=std::max(cl(0,k),
                                        std::min(cl(1,k),
                                            std::copysign(au,u)/(v(k)+vp(k)*al2)));
                            }
                            if (a(k) == ak) continue; 
                            auto d=a(k)-ak; 
                            dlx=std::max(dlx,v(k)*d*d);

                            // TODO
                            auto d_scaled = d/xs(k);
                            wr -= d_scaled*x.col(k).cwiseProduct(w);
                            uu -= d_scaled * xb(k);
                            tt -= d_scaled * xm(k);
                            //wr.array()-=d*w.array()*x.col(k).array(); 
                            //f += d*x.col(k);
                        }
                        if (intr) { 
                            // TODO
                            //auto d=wr.sum()/v0; 
                            auto d = tt/v0 - uu;

                            az += d;
                            dlx=std::max(dlx,v0*d*d); 

                            // TODO
                            uu += d;
                            //wr-=d*w; 
                            //f.array()+=d;
                        }
                        if (dlx < thr) break; 
                        if (nlp > maxit) throw util::maxit_reached_error();
                    }
                }

                // TODO
                t.array() = g.array();
                for (int i = 0; i < nin; ++i) {
                    t += (a(ia(i)-1) / xs(ia(i)-1)) * x.col(ia(i)-1);
                }
                if (apply_suggested_fixes) {
                    w.array() = q.array() *
                        ((t.array()+uu).abs().min(fmax)).matrix().binaryExpr(t,
                            [&](auto x, auto y) { return std::copysign(x,y+uu); }).array().exp();
                } else {
                    auto euu = std::exp(std::copysign(std::min(std::abs(uu), fmax), uu));
                    //w.array() = q.array() * 
                    //    (f.array().abs().min(fmax)).binaryExpr(f.array(),
                    //        [](auto x, auto y) { return std::copysign(x,y); }).array().exp();
                    w.array() = euu * q.array() *
                        (t.array().abs().min(fmax)).matrix().binaryExpr(t,
                            [](auto x, auto y) { return std::copysign(x,y); }).array().exp();
                }
                v0 = w.sum(); 
                //wr=t-w;
                wr = qy - w*(static_cast<FloatType>(1.0)-uu);
                tt = wr.sum(); 

                auto d = az - az0;
                if (v0*d*d < thr) {
                    bool ix = false;
                    for (int j = 0; j < nin; ++j) {
                        auto k = ia(j)-1;
                        auto dd = a(k) - as(k);
                        if (v(k) * dd * dd < thr) continue;
                        ix = true;
                        break;
                    }
                    if (!ix) {
                        for (int k = 0; k < ni; ++k) {
                            if (ixx(k) || !ju[k]) continue;

                            // TODO
                            xm(k) = x.col(k).dot(w);
                            auto gk = (x.col(k).dot(wr)-uu*(xm(k)-v0*xb(k))-xb(k)*tt)/xs(k);
                            ga(k) = std::abs(gk);

                            if (ga(k) > al1 * vp(k)) {
                                ixx[k] = true;
                                ix = true;
                            }
                        }
                        if (!ix) break;
                    }
                }
            } 
        } 
        catch (const util::maxit_reached_error& e) {
            jerr = e.err_code(ilm);
            return;
        }
        catch (const util::elnet_error& e) {
            jerr = e.err_code(ilm);
            break;
        }

        for (int j = 0; j < nin; ++j) {
            ca(j, ilm) = a(ia(j)-1);
        }
        kin(ilm) = nin;
        a0(ilm) = az;
        alm(ilm) = al;
        lmu=ilm+1;

        // TODO
        //dev(ilm) = (t.dot(f) - v0 - dv0) / dvr;
        dev(ilm) = (qy.dot(t) + yb*uu -v0 -dv0)/dvr;

        if (lmu < mnl) continue;
        if (flmin >= 1.0) continue;
        auto me = (ca.col(ilm).head(nin).array() != 0).count();
        if (me > ne) break;
        if ((dev(ilm) - dev(ilm-mnl+1)) < sml * dev(ilm)) break;
        if (dev(ilm) > devmax) break;
    }

    // TODO
    g.array() = t.array() + uu;
}

} // namespace transl
} // namespace glmnetpp
