#pragma once
#include <limits>
#include <testutil/mock_pb.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>
#include <testutil/internal.hpp>
#include <Eigen/Core>

namespace glmnetpp {
namespace transl {

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
        , class CAType
        , class IAType
        , class KinType
        , class DevType
        , class AlmType>
inline void fishnet1(
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
    Eigen::VectorXd f(no);

    auto bta = beta; 
    auto omb = 1.0 - bta;
    t.array()=q.array()*y.array(); 
    auto yb=t.sum(); 
    auto fmax=std::log(std::numeric_limits<double>::max()*0.1);
    double dv0 = 0;
    double az = 0;
    if ((g.array() == 0).all()) {
        if (intr) { 
            w=q*yb;  
            az=std::log(yb); 
            dv0=yb*(az-1.0); 
            f.array() = az;
        }
        else { 
            w=q; 
            az = 0;
            dv0=-1.0;
            f.setZero();
        }
    }
    else {
        w.array() = q.array() * (g.array().abs().min(fmax)).matrix().binaryExpr(g,
                    [&](auto x, auto y) { return std::copysign(x, y); }).array().exp();
        auto v0=w.sum();
        f = g;
        if (intr) { 
            auto eaz=yb/v0; 
            w *= eaz; 
            az=std::log(eaz); 
            dv0=t.dot(g)-yb*(1.0-az);
            f.array() += az;
        }
        else { 
            az=0.0; 
            dv0=t.dot(g)-v0; }
    }

    wr=t-w; 
    auto v0=1.0; 
    if (intr) v0=yb; 
    auto dvr=-yb;
    for (int i = 0; i < no; ++i) {
        if (t(i) > 0.0) dvr += t(i)*std::log(y(i));
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
        ga(j)=std::abs(wr.dot(x.col(j)));
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
        auto tlam = bta*(2.0*al-al0);
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
                    if (ixx(j)) v(j) = w.dot(x.col(j).array().square().matrix());
                }

                // CD loop
                while (1) {
                    ++nlp; 
                    auto dlx=0.0;
                    for (int k = 0; k < ni; ++k) {
                        if (!ixx(k)) continue; 
                        auto ak=a(k);
                        auto gk = wr.dot(x.col(k));
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
                        wr.array()-=d*w.array()*x.col(k).array(); 
                        if (!mm(k)) { 
                            ++nin; 
                            if (nin > nx) throw util::max_active_reached_error();
                            mm(k)=nin; 
                            ia(nin-1)=k+1;
                        }
                    }
                    if (intr) { 
                        auto d=wr.sum()/v0;
                        az += d; 
                        dlx = std::max(dlx,v0*d*d); 
                        wr-=d*w; 
                    }
                    if (dlx < thr) break; 
                    if (nlp > maxit) throw util::maxit_reached_error(); 
                    while (1) { 
                        ++nlp; 
                        auto dlx=0.0;
                        for (int l = 0; l < nin; ++l) {
                            auto k=ia(l)-1; 
                            auto ak=a(k);
                            auto gk = wr.dot(x.col(k));
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
                            wr.array()-=d*w.array()*x.col(k).array(); 
                        }
                        if (intr) { 
                            auto d=wr.sum()/v0; 
                            az += d;
                            dlx=std::max(dlx,v0*d*d); 
                            wr-=d*w; 
                        }
                        if (dlx < thr) break; 
                        if (nlp > maxit) throw util::maxit_reached_error();
                    }
                }

                f.array() = az + g.array();
                for (int i = 0; i < nin; ++i) {
                    auto k = ia(i)-1;
                    f += a(k) * x.col(k);
                }
                w.array() = q.array() * 
                    (f.array().abs().min(fmax)).binaryExpr(f.array(),
                        [](auto x, auto y) { return std::copysign(x,y); }).array().exp();
                v0 = w.sum(); 
                wr=t-w;
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
                            ga(k) = std::abs(wr.dot(x.col(k)));
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
        dev(ilm) = (t.dot(f) - v0 - dv0) / dvr;
        if (lmu < mnl) continue;
        if (flmin >= 1.0) continue;
        auto me = (ca.col(ilm).head(nin).array() != 0).count();
        if (me > ne) break;
        if ((dev(ilm) - dev(ilm-mnl+1)) < sml * dev(ilm)) break;
        if (dev(ilm) > devmax) break;
    }

    g = f;
}

} // namespace transl
} // namespace glmnetpp
