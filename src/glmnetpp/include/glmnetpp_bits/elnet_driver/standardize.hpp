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

struct MultStandardize1
{
    template <class XType
            , class YType
            , class WType
            , class JUType
            , class XMType
            , class XSType
            , class YMType
            , class YSType
            , class XVType
            , class ValueType
            >
    static void eval(
        XType& x, 
        YType& y, 
        WType& w, 
        bool isd, 
        bool jsd,
        bool intr, 
        const JUType& ju, 
        XMType& xm,
        XSType& xs,
        YMType& ym,
        YSType& ys,
        XVType& xv,
        ValueType& ys0)
    {
        using vec_t = Eigen::Matrix<ValueType, Eigen::Dynamic, 1>;

        auto ni = x.cols();
        auto nr = y.cols();

        w /= w.sum();
        vec_t v = w.array().sqrt().matrix();

        if (!intr) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue;
                xm(j) = 0.0;
                x.col(j).array() *= v.array();
                auto z = x.col(j).squaredNorm();
                if (isd) {
                    auto xbq = v.dot(x.col(j));
                    xbq *= xbq;
                    auto vc = z-xbq;
                    xs(j) = std::sqrt(vc);
                    x.col(j) /= xs(j);
                    xv(j) = 1.0 + xbq / vc;
                }
                else {
                    xs(j) = 1.0;
                    xv(j) = z;
                }
            }
            ys0 = 0.0;
            for (int j = 0; j < nr; ++j) {
                ym(j) = 0.0;
                y.col(j).array() *= v.array();
                auto z = y.col(j).squaredNorm();
                if (jsd) {
                    auto u = v.dot(y.col(j));
                    u = z - u*u;
                    ys0 += z/u;
                    ys(j) = std::sqrt(u);
                    y.col(j) /= ys(j);
                }
                else {
                    ys(j) = 1.0;
                    ys0 += z;
                }
            }
            return;
        }

        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = w.dot(x.col(j));
            x.col(j).array() = v.array() * (x.col(j).array()-xm(j));
            xv(j) = x.col(j).squaredNorm();
            if (isd) xs(j) = std::sqrt(xv(j));
        }

        if (!isd) { xs.setOnes(); }
        else {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue;
                x.col(j) /= xs(j);
            }
            xv.setOnes();
        }

        ys0 = 0.0;
        for (int j = 0; j < nr; ++j) {
            ym(j) = w.dot(y.col(j));
            y.col(j).array() = v.array() * (y.col(j).array() - ym(j));
            auto z = y.col(j).squaredNorm();
            if (jsd) { 
                ys(j) = std::sqrt(z);  
                y.col(j) /= ys(j);
            } else {
                ys0 += z;
            }
        }
        if (!jsd) ys.setOnes();
        else ys0 = nr;
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

struct MultSpStandardize1
{
    template <class XType
            , class YType
            , class WType
            , class JUType
            , class XMType
            , class XSType
            , class YMType
            , class YSType
            , class XVType
            , class ValueType
            >
    static void eval(
        const XType& x, 
        YType& y, 
        WType& w, 
        bool isd, 
        bool jsd,
        bool intr, 
        const JUType& ju, 
        XMType& xm,
        XSType& xs,
        YMType& ym,
        YSType& ys,
        XVType& xv,
        ValueType& ys0)
    {
        auto ni = x.cols();
        auto nr = y.cols();

        w /= w.sum();

        if (!intr) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue;
                xm(j) = 0.0;
                auto xsq = x.col(j).cwiseProduct(x.col(j));
                auto z = xsq.dot(w);
                if (isd) {
                    auto xbq = x.col(j).dot(w);
                    xbq *= xbq;
                    auto vc = z-xbq;
                    xs(j) = std::sqrt(vc);
                    xv(j) = 1.0 + xbq / vc;
                }
                else {
                    xs(j) = 1.0;
                    xv(j) = z;
                }
            }
            ys0 = 0.0;
            for (int j = 0; j < nr; ++j) {
                ym(j) = 0.0;
                auto ysq = y.col(j).cwiseProduct(y.col(j));
                auto z = ysq.dot(w);
                if (jsd) {
                    auto u = y.col(j).dot(w);
                    u = z - u*u;
                    ys0 += z/u;
                    ys(j) = std::sqrt(u);
                    y.col(j) /= ys(j);
                }
                else {
                    ys(j) = 1.0;
                    ys0 += z;
                }
            }
            return;
        }

        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = x.col(j).dot(w);
            xv(j) = x.col(j).cwiseProduct(x.col(j)).dot(w) - xm(j)*xm(j);
            if (isd) xs(j) = std::sqrt(xv(j));
        }

        if (!isd) { xs.setOnes(); }
        else {
            xv.setOnes();
        }

        ys0 = 0.0;
        for (int j = 0; j < nr; ++j) {
            ym(j) = w.dot(y.col(j));
            y.col(j).array() -= ym(j);
            auto z = y.col(j).array().square().matrix().dot(w);
            if (jsd) { 
                ys(j) = std::sqrt(z);  
                y.col(j) /= ys(j);
            } else {
                ys0 += z;
            }
        }
        if (!jsd) ys.setOnes();
        else ys0 = nr;
    }
};

struct LStandardize1
{

    template <class XType
            , class WType
            , class JUType
            , class IntType
            , class XMType
            , class XSType>
    static void eval(
            XType& x,
            const WType& w,
            const JUType& ju,
            IntType isd,
            IntType intr,
            XMType& xm,
            XSType& xs
            )
    {
        auto ni = x.cols();

        // without intercept
        if (!intr) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue; 
                xm(j) = 0.0;
                if (isd) { 
                    auto mean = w.dot(x.col(j));
                    auto vc = w.dot(x.col(j).array().square().matrix())
                        - mean * mean;
                    xs(j) = std::sqrt(vc); 
                    x.col(j) /= xs(j);
                }
            }
            return;
        }

        // with intercept
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = w.dot(x.col(j)); 
            x.col(j).array() -= xm(j);
            if (isd) { 
                xs(j) = std::sqrt(w.dot(x.col(j).array().square().matrix())); 
                x.col(j) /= xs(j);
            }
        }
    }
};

struct MultLStandardize1
{

    template <class XType
            , class WType
            , class JUType
            , class IntType
            , class XMType
            , class XSType
            , class XVType>
    static void eval(
            XType& x,
            const WType& w,
            const JUType& ju,
            IntType isd,
            IntType intr,
            XMType& xm,
            XSType& xs,
            XVType& xv
            )
    {
        auto ni = x.cols();

        // without intercept
        if (!intr) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue; 
                xm(j) = 0.0;
                xv(j) = x.col(j).array().square().matrix().dot(w);
                if (isd) { 
                    auto mean_sq = w.dot(x.col(j));
                    mean_sq *= mean_sq; // now it is truly mean-squared
                    auto vc = xv(j) - mean_sq;
                    xs(j) = std::sqrt(vc); 
                    x.col(j) /= xs(j);
                    xv(j) = 1.0 + mean_sq/vc;
                }
            }
            return;
        }

        // with intercept
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = w.dot(x.col(j)); 
            x.col(j).array() -= xm(j);
            xv(j) = x.col(j).array().square().matrix().dot(w);
            if (isd) { 
                xs(j) = std::sqrt(xv(j));
                x.col(j) /= xs(j);
                xv(j) = 1.0;
            }
        }
    }
};

struct SpLStandardize2
{

    template <class XType
            , class WType
            , class JUType
            , class IntType
            , class XMType
            , class XSType>
    static void eval(
            const XType& x,
            const WType& w,
            const JUType& ju,
            IntType isd,
            IntType intr,
            XMType& xm,
            XSType& xs
            )
    {
        auto ni = x.cols();

        // without intercept
        if (!intr) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue; 
                xm(j) = 0.0;
                if (isd) { 
                    auto vc = x.col(j).cwiseProduct(x.col(j)).dot(w);
                    auto mean = x.col(j).dot(w);
                    vc -= mean * mean;
                    xs(j) = std::sqrt(vc); 
                }
                else { xs(j) = 1.0; }
            }
            return;
        }

        // with intercept
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = x.col(j).dot(w); 
            if (isd) { 
                xs(j) = std::sqrt(
                        x.col(j).cwiseProduct(x.col(j)).dot(w)
                         - xm(j) * xm(j)); 
            }
        }
        if (!isd) xs.array() = 1.0;
    }
};

struct MultSpLStandardize2
{
    template <class XType
            , class WType
            , class JUType
            , class IntType
            , class XMType
            , class XSType
            , class XVType>
    static void eval(
            const XType& x,
            const WType& w,
            const JUType& ju,
            IntType isd,
            IntType intr,
            XMType& xm,
            XSType& xs,
            XVType& xv
            )
    {
        auto ni = x.cols();

        // without intercept
        if (!intr) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) continue; 
                xm(j) = 0.0;
                xv(j) = x.col(j).cwiseProduct(x.col(j)).dot(w);
                if (isd) { 
                    auto mean_sq = x.col(j).dot(w);
                    mean_sq *= mean_sq; // now mean squared
                    auto vc = xv(j) - mean_sq;
                    xs(j) = std::sqrt(vc); 
                    xv(j) = 1.0 + mean_sq/vc;
                }
                else { xs(j) = 1.0; }
            }
            return;
        }

        // with intercept
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) continue;
            xm(j) = x.col(j).dot(w); 
            xv(j) = x.col(j).cwiseProduct(x.col(j)).dot(w) - xm(j)*xm(j);
            if (isd) { 
                xs(j) = std::sqrt(xv(j)); 
                xv(j) = 1.0;
            }
        }
        if (!isd) xs.array() = 1.0;
    }
};

} // namespace glmnetpp
