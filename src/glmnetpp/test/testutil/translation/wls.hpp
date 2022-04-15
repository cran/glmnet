#include <cstddef>
#include <Eigen/Core>
#include <vector>
#include <cmath>
#include <glmnetpp_bits/util/type_traits.hpp> 
#include <glmnetpp_bits/util/iterator/counting_iterator.hpp>
#include <glmnetpp_bits/util/exceptions.hpp>

namespace glmnetpp {
namespace transl {
namespace details {

// Non-dense version.
template <class XType, class=void>
struct Square
{
    static auto eval(const XType& x) 
    {
        return x.cwiseProduct(x);
    }
};

// Dense version.
template <class XType>
struct Square<XType, 
    typename std::enable_if<util::is_dense<XType>::value>::type>
{
    static auto eval(const XType& x)
    {
        return x.array().square().matrix();
    }
};

// Some helper utilities that need to be distinguished 
// between sparse and dense matrices.
template <class XType>
inline auto square(const XType& x)
{ return Square<XType>::eval(x); }

template <class XType, class RType, class ValueType, class=void>
struct Dot
{
    static auto eval(const XType& x, const RType& r, ValueType r_sum)
    {
        auto d = x.dot(r);
        d -= r_sum * x.mean();
        d /= x.sd(); 
        return d;
    }
};

template <class XType, class RType, class ValueType>
struct Dot<XType, RType, ValueType, 
    typename std::enable_if<util::is_dense<XType>::value>::type >
{
    static auto eval(const XType& x, const RType& r, ValueType)
    {
        return x.dot(r);
    }
};

// Compute the dot product between r and x vectors.
// If x is sparse, it must support dot(r, r_sum),
// where r_sum is the sum of all entries of r.
//
// @param   r_sum   sum of all entries in r. It is more efficient to supply 
//                  this argument to compute the dot product for sparse x. 
//                  It is ignored if x is dense.
template <class XType, class RType, class ValueType>
inline auto dot(
    const XType& x,
    const RType& r,
    ValueType r_sum)
{ return Dot<XType, RType, ValueType>::eval(x, r, r_sum); }

template <class XType, class VType, class ValueType, class=void>
struct ComputeXV
{
    static auto eval(const XType& x, const VType& v, ValueType v_sum)
    {
        auto xv = square(x).dot(v);
        xv -= 2 * x.mean() * x.dot(v);
        return (xv + v_sum * x.mean() * x.mean()) / (x.sd() * x.sd()); 
    }
};

template <class XType, class VType, class ValueType>
struct ComputeXV<XType, VType, ValueType,
    typename std::enable_if<util::is_dense<XType>::value>::type >
{
    static auto eval(const XType& x, const VType& v, ValueType)
    {
        return square(x).dot(v);
    }
};

// Compute xv temporary helper.
template <class XType, class VType, class ValueType>
inline auto compute_xv(const XType& x, 
                       const VType& v,
                       ValueType v_sum)
{ return ComputeXV<XType, VType, ValueType>::eval(x,v,v_sum); }

template <class RType, class XType, class VType, class ValueType, class=void>
struct UpdateR
{
    static auto eval(RType& r, const XType& x, const VType& v, ValueType d)
    {
        auto d_scaled = d / x.sd();
        r -= d_scaled * x.cwiseProduct(v);
        r += (d_scaled * x.mean()) * v;
    }
};

template <class RType, class XType, class VType, class ValueType>
struct UpdateR<RType, XType, VType, ValueType,
    typename std::enable_if<util::is_dense<XType>::value>::type >
{
    static auto eval(RType& r, const XType& x, const VType& v, ValueType d)
    {
        r -= (d * v.array() * x.array()).matrix();
    }
};


// Update residual vector
template <class RType
        , class XType
        , class VType
        , class ValueType>
inline void update_r(RType& r, 
                     const XType& x,
                     const VType& v,
                     ValueType d)
{ 
    using r_t = typename std::decay<RType>::type;
    UpdateR<r_t, XType, VType, ValueType>::eval(r, x, v, d); 
}

template <bool do_update, class ValueType, class RType, class=void>
struct UpdateRSum
{
    static_assert(!do_update, "Not a valid parameter list for updating r_sum.");
    static void eval(ValueType&, const RType&) {}
};

template <class ValueType, class RType>
struct UpdateRSum<false, ValueType, RType, void>
{
    static void eval(ValueType&, const RType&) {}
};

template <class ValueType, class RType>
struct UpdateRSum<true, ValueType, RType,
    typename std::enable_if<util::is_dense<RType>::value>::type >
{
    static void eval(ValueType& r_sum, const RType& r)
    {
        r_sum = r.sum();
    }
};

template <class ValueType, class RType>
struct UpdateRSum<true, ValueType, RType,
    typename std::enable_if<std::is_convertible<RType, ValueType>::value>::type >
{
    static void eval(ValueType& r_sum, const RType& r)
    {
        r_sum = r;
    }
};

template <bool do_update
        , class ValueType
        , class RType>
inline void update_r_sum(ValueType& r_sum, const RType& r)
{ 
    using value_t = typename std::decay<ValueType>::type;
    UpdateRSum<do_update, value_t, RType>::eval(r_sum, r);
}

template <bool partial_update>
struct GetUpdateRSumVal
{
    template <class ValueType, class RType>
    static auto eval(ValueType diff, const RType&)
    {
        return diff;
    }
};

template <>
struct GetUpdateRSumVal<false>
{
    template <class ValueType, class RType>
    static auto eval(ValueType, const RType& r)
    {
        return r;
    }
};

// Update intercept term (and other invariant quantities)
template <bool do_update_r_sum
        , bool partial_update
        , class IntType
        , class ValueType
        , class VType
        , class RType>
inline void update_intercept(
        IntType intr,
        ValueType xmz,
        const VType& v,
        RType& r,
        ValueType& r_sum,
        ValueType& aint,
        ValueType& rsqc,
        ValueType& dlx)
{
    double d = 0.0; 
    if (intr != 0) { 
        // if we haven't been updating r_sum,
        // we must at least update here for the next if-block.
        update_r_sum<!do_update_r_sum>(r_sum, r);
        d = r_sum / xmz; 
    }
    if (d != 0.0) {
        aint += d;
        auto diff = r_sum - d * xmz;
        rsqc += d * (r_sum + diff);
        dlx = std::max(dlx, xmz * d * d);
        r -= d * v;
        auto&& update_val = GetUpdateRSumVal<partial_update>::eval(diff, r);
        update_r_sum<do_update_r_sum>(r_sum, update_val);
    }
}

template <bool do_active>
struct AddActive
{
    template <class IntType, class AddActiveFType>
    static void eval(IntType k, AddActiveFType add_active_f)
    {
        add_active_f(k);
    }
};

template <>
struct AddActive<false>
{
    template <class IntType, class AddActiveFType>
    static void eval(IntType, AddActiveFType) {}
};

template <bool do_update_r_sum
        , bool add_active
        , class Iter
        , class XType
        , class VType
        , class RType
        , class ValueType
        , class AType
        , class XVType
        , class VPType
        , class CLType
        , class SkipFType
        , class AddActiveFType>
inline void coord_desc(
        Iter begin,
        Iter end,
        const XType& x,
        const VType& v,
        RType& r,
        ValueType& r_sum,
        AType& a,
        const XVType& xv,
        const VPType& vp,
        ValueType ab,
        ValueType dem,
        const CLType& cl,
        ValueType& rsqc,
        ValueType& dlx,
        SkipFType skip_f,
        AddActiveFType add_active_f)
{
    for (; begin != end; ++begin) {
        
        auto k = *begin;
        
        if (skip_f(k)) continue;
        
        // check if ST threshold for descent is met
        // if it goes to 0, set a(k) = 0.0
        // if not, set a(k) to the post-gradient descent value
        // u is the kth partial residual
        auto gk = details::dot(x.col(k), r, r_sum);
        auto ak = a(k); 
        auto u = gk + ak * xv(k); 
        auto au = std::abs(u) - vp(k) * ab;
        if (au < 0.0) { a(k) = 0.0; }
        else {
            a(k) = std::max(cl(0,k),
                    std::min(cl(1,k), 
                        std::copysign(au, u) / (xv(k) + vp(k) * dem)
                        ) );
        }
        
        // if the update didn't change the coefficient value, go to
        // the next feature/variable
        if (a(k) == ak) continue;
        
        AddActive<add_active>::eval(k, add_active_f);

        // update residual r, rsqc, and dlx (diff exp from wls)
        auto d = a(k) - ak;
        rsqc += d * (2.0 * gk - d * xv(k));
        details::update_r(r, x.col(k), v, d);
        details::update_r_sum<do_update_r_sum>(r_sum, r);
        dlx = std::max(xv(k) * d * d, dlx);
            
    }
}

// Partially fit WLS on the active set.
template <bool do_update_r_sum
        , class IntType
        , class ValueType
        , class IAType
        , class RType
        , class XMapType
        , class VType
        , class AType
        , class XVType
        , class VPType
        , class CLType>
inline void wls_partial_fit(
        IntType& iz, 
        IntType& jz,
        IntType& nlp,
        IntType& nino,
        IntType maxit,
        IntType intr,
        ValueType thr,
        ValueType xmz,
        IAType& ia,
        RType& r,
        ValueType& r_sum,
        const XMapType& x,
        const VType& v,
        AType& a,
        ValueType& aint,
        const XVType& xv,
        const VPType& vp,
        const CLType& cl,
        ValueType& rsqc,
        ValueType ab,
        ValueType dem)
{
    using value_t = ValueType;

    iz = 1;

    while (1) {
        
        ++nlp;
        value_t dlx = 0.0;
        coord_desc<do_update_r_sum, false>(
            ia.data(), ia.data() + nino,
            x, v, r, r_sum, a, xv, vp, ab, dem, cl, rsqc, dlx,
            [](auto) { return false; },
            [](auto) {}
        );

        // updating of intercept term   
        update_intercept<do_update_r_sum, true>(
                intr, xmz, v, r, r_sum, aint, rsqc, dlx);

        if (dlx < thr) break;
        if (nlp > maxit) { 
            throw util::maxit_reached_error();
        }
    }

    // set jz = 0 so that we have to go into :again: tag
    // to check KKT conditions.
    jz = 0;
}

} // namespace details

/*
 * Experimental C++ implementation of WLS
 *
 * alm0: previous lambda value
 * almc: current lambda value
 * alpha: alpha
 * m: current lambda iteration no.
 * no: no of observations
 * ni: no of variables
 * x: x matrix
 * r: weighted residual! v * (y - yhat)
 * v: weights
 * intr: 1 if fitting intercept, 0 otherwise
 * ju: ju(k) = 1 if feature k is included for consideration in model
 * vp: relative penalties on features (sum to ni)
 * cl: coordinate limits
 * nx: limit for max no. of variables ever to be nonzero
 * thr: threshold for dlx
 * maxit: max no. of passes over the data for all lambda values
 * a, aint: warm start for coefs, intercept
 * g: abs(dot_product(r,x(:,j)))
 * ia: mapping nino to k
 * iy: ever-active set (for compatibility with sparse version)
 * iz: flag for loop. 0 for first lambda, 1 subsequently
 * mm: mapping k to nino
 * nino: no. of features that have ever been nonzero
 * rsqc: R^2
 * nlp: no. of passes over the data
 * jerr: error code
 */

template <class ValueType
        , class IntType
        , class XMapType
        , class RType
        , class XVType
        , class VType
        , class JUType
        , class VPType
        , class CLType
        , class AType
        , class GType
        , class IAType
        , class IYType
        , class MMType>
inline void wls(
    ValueType alm0,
    ValueType almc,
    ValueType alpha,
    IntType m,
    const XMapType& x,
    RType& r,
    XVType& xv,
    const VType& v,
    IntType intr,
    const JUType& ju,
    const VPType& vp,
    const CLType& cl,
    IntType nx,
    ValueType thr,
    IntType maxit,
    AType& a,
    ValueType& aint,
    GType& g,
    IAType& ia,
    IYType& iy,
    IntType& iz,
    MMType& mm,
    IntType& nino,
    ValueType& rsqc,
    IntType& nlp,
    IntType& jerr)
{
    using value_t = ValueType;
    using int_t = IntType;

    auto ni = x.cols();

    // compute xmz (for intercept later)
    const value_t xmz = v.sum();
    
    // only update r_sum if X is not dense.
    constexpr bool do_update_r_sum = !util::is_dense<XMapType>::value;
    value_t r_sum = 0.0;
    details::update_r_sum<do_update_r_sum>(r_sum, r);
    
    // compute g initialization
    for (int_t j = 0; j < ni; ++j) {
        if (ju(j) == 0) continue; 
        g(j) = std::abs(details::dot(x.col(j), r, r_sum));
    }

    // compute xv
    // Note: weight v may have changed from previous call to wls,
    // so we must update existing xv values.
    for (int_t j = 0; j < ni; ++j) {
        if (iy(j) > 0) {
            xv(j) = details::compute_xv(x.col(j), v, xmz);
        }
    }

    // ab: lambda * alpha, dem: lambda * (1 - alpha)
    value_t ab = almc * alpha; 
    value_t dem = almc * (1.0 - alpha);

    // strong rules: iy(k) = 1 if we don't discard feature k
    value_t tlam = alpha * (2.0 * almc - alm0);
    for (int_t k = 0; k < ni; ++k) {
        if (iy(k) == 1 || ju(k) == 0) continue; 
        if (g(k) > tlam * vp(k)) {
            iy(k) = 1; 
            xv(k) = details::compute_xv(x.col(k), v, xmz);
        }
    }
    
    int_t jz = 1;

    try {
        if (iz*jz != 0) {
            details::wls_partial_fit<do_update_r_sum>(
                iz, jz, nlp, nino, maxit, intr,
                thr, xmz, ia, r, r_sum, x, v, a, aint, xv, 
                vp, cl, rsqc, ab, dem);
        }
        
        while (1) {

            bool converged = false;
            while (1) {

                // :again: 
                ++nlp; 
                value_t dlx = 0.0;

                details::coord_desc<do_update_r_sum, true>(
                    util::counting_iterator<int_t>(0),
                    util::counting_iterator<int_t>(ni),
                    x, v, r, r_sum, a, xv, vp, ab, dem, cl, rsqc, dlx,
                    [&](auto k) { return iy(k) == 0; },
                    [&, nx](auto k) {
                        // if coef for feature k was previously 0, we now have a 
                        // new non-zero coef. update nino, mm(k) and ia(nino).
                        // This is only needed if we are _not_ iterating over active set
                        // because mm(k) will always be > 0 if k is in active set, by def.
                        if (mm(k) == 0) {
                            ++nino;
                            if (nino > nx) 
                                throw util::max_active_reached_error();
                            mm(k) = nino; 
                            // Note: ia is never used outside this function call,
                            // so it is OK to store indices in 0-index.
                            ia(nino-1) = k+1; 
                        }
                    }
                );

                // updating of intercept term
                details::update_intercept<do_update_r_sum, false>(
                        intr, xmz, v, r, r_sum, aint, rsqc, dlx);
                
                // in wls, this leads to KKT checks. here, we exit
                // the loop instead.
                if (dlx < thr) {
                    bool ixx = false;
                    for (int_t k = 0; k < ni; ++k) {
                        if (iy(k) == 1 || ju(k) == 0) continue; 
                        g(k) = std::abs(details::dot(x.col(k), r, r_sum));
                        if (g(k) > ab * vp(k)) {
                            iy(k) = 1; 
                            xv(k) = details::compute_xv(x.col(k), v, xmz);
                            ixx = true;
                        }
                    }
                    if (!ixx) {
                        converged = true;
                        break;
                    }
                }

                else break;

            } // end :again: while

            if (converged) break;

            // if we've gone over max iterations, return w error
            if (nlp > maxit) { 
                throw util::maxit_reached_error();
            }

            // this is like the :b: loop in wls (M)
            details::wls_partial_fit<do_update_r_sum>(
                iz, jz, nlp, nino, maxit, intr,
                thr, xmz, ia, r, r_sum, x, v, a, aint, xv, 
                vp, cl, rsqc, ab, dem);

        } // end outer while
    }
    catch (const util::elnet_error& e) {
        jerr = e.err_code(m-1); // m is 1-indexed; err_code expects 0-indexed
    }
}

} // namespace transl
} // namespace glmnetpp
