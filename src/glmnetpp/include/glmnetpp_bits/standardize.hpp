#pragma once
#include <Eigen/Core>
#include <type_traits>

namespace glmnetpp {

struct Standardize1
{
    template <class XType
            , class YType
            , class WType
            , class JUType
            , class XMType
            , class XSType
            , class ValueType
            , class XVType>
    static void eval(
        XType& x, 
        YType& y, 
        WType& w, 
        bool isd, 
        bool intr, 
        const JUType& ju, 
        XMType& xm,
        XSType& xs,
        ValueType& ym,
        ValueType& ys,
        XVType& xv)
    {
        using value_t = typename std::decay_t<XType>::Scalar;
        using vec_t = typename Eigen::Matrix<value_t, Eigen::Dynamic, 1>;

        auto ni = x.cols();

        w /= w.sum(); 

        vec_t v = w.array().sqrt().matrix();

        // without intercept
        if (!intr) {

            ym = 0.0;
            y.array() *= v.array();

            // trevor changed 3/24/2020
            ys = y.norm(); 
            y /= ys;
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue; 
                xm(j) = 0.0; 
                auto x_j = x.col(j);
                x_j.array() *= v.array();
                xv(j) = x_j.squaredNorm();
                if (isd) { 
                    auto xbq = x_j.dot(v);
                    xbq *= xbq;
                    auto vc = xv(j) - xbq;
                    xs(j) = std::sqrt(vc); 
                    x_j /= xs(j); 
                    xv(j) = 1.0 + xbq/vc;
                }
                else { xs(j) = 1.0; }
            }

        } 

        // with intercept
        else {

            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue;
                auto x_j = x.col(j);
                xm(j) = x_j.dot(w); 
                x_j.array() = v.array() * (x_j.array() - xm(j));
                xv(j) = x_j.squaredNorm(); 
                if (isd) xs(j) = std::sqrt(xv(j));
            }
            if (!isd) { xs.setOnes(); }
            else {
                 for (int j = 0; j < ni; ++j) {
                     if (!ju[j]) continue; 
                     auto x_j = x.col(j);
                     x_j /= xs(j);
                 }
                 xv.setOnes();
            }
            ym = y.dot(w);
            y.array() = v.array() * (y.array() - ym); 
            ys = y.norm(); 
            y /= ys;
        }
    }
};

struct Standardize
{

    template <class XType
            , class YType
            , class WType
            , class JUType
            , class GType
            , class XMType
            , class XSType
            , class ValueType
            , class XVType>
    static void eval(
        XType& x, 
        YType& y, 
        WType& w, 
        bool isd, 
        bool intr, 
        const JUType& ju, 
        GType& g,
        XMType& xm,
        XSType& xs,
        ValueType& ym,
        ValueType& ys,
        XVType& xv)
    {
        auto ni = x.cols();
        Standardize1::eval(x, y, w, isd, intr, ju, xm, xs, ym, ys, xv);
        g.setZero(); 
        for (int j = 0; j < ni; ++j) {
            if (ju[j]) g(j) = x.col(j).dot(y);
        }
    }
};

struct SpStandardize1
{
    template <class XType
            , class YType
            , class WType
            , class JUType
            , class XMType
            , class XSType
            , class ValueType
            , class XVType>
    static void eval(
        const XType& x, 
        YType& y, 
        WType& w, 
        bool isd, 
        bool intr, 
        const JUType& ju, 
        XMType& xm,
        XSType& xs,
        ValueType& ym,
        ValueType& ys,
        XVType& xv)
    {
        auto ni = x.cols();
        w /= w.sum();

        // without intercept
        if (!intr) {
            ym = 0.0;

            // trevor changed 3/24/2020
            ys = std::sqrt(y.array().square().matrix().dot(w)); 
            y /= ys;
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue; 
                xm(j) = 0.0; 
                auto x_j = x.col(j);
                xv(j) = x_j.cwiseProduct(x_j).dot(w);
                if (isd) { 
                    auto xbq = x_j.dot(w);
                    xbq *= xbq;
                    auto vc = xv(j) - xbq;
                    xs(j) = std::sqrt(vc); 
                    xv(j) = 1.0 + xbq/vc;
                }
                else { xs(j) = 1.0; }
            }
        }

        // with intercept
        else {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue;
                auto x_j = x.col(j);
                xm(j) = x_j.dot(w); 
                xv(j) = x_j.cwiseProduct(x_j).dot(w) - xm(j) * xm(j);
                if (isd) xs(j) = std::sqrt(xv(j));
            }
            if (!isd) { xs.setOnes(); }
            else { xv.setOnes(); }
            ym = y.dot(w);
            y.array() -= ym; 
            ys = std::sqrt(y.array().square().matrix().dot(w)); 
            y /= ys;
        }
    }
};

struct SpStandardize
{

    template <class XType
            , class YType
            , class WType
            , class JUType
            , class GType
            , class XMType
            , class XSType
            , class ValueType
            , class XVType>
    static void eval(
        const XType& x, 
        YType& y, 
        WType& w, 
        bool isd, 
        bool intr, 
        const JUType& ju, 
        GType& g,
        XMType& xm,
        XSType& xs,
        ValueType& ym,
        ValueType& ys,
        XVType& xv)
    {
        auto ni = x.cols();
        SpStandardize1::eval(x, y, w, isd, intr, ju, xm, xs, ym, ys, xv);
        g.setZero();
        for (int j = 0; j < ni; ++j) {
            if (ju[j]) {
                g(j) = x.col(j).dot(
                        (y.array() * w.array()).matrix() ) / xs(j);
            }
        }
    }
};

} // namespace glmnetpp
