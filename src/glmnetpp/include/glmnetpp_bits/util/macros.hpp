#pragma once

/*
 * GLMNETPP_STRONG_INLINE is a stronger version of the inline, 
 * using __forceinline on MSVC, always_inline on GCC/clang, and otherwise just use inline.
 */
#ifndef GLMNETPP_STRONG_INLINE
#if defined(_MSC_VER)
#define GLMNETPP_STRONG_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define GLMNETPP_STRONG_INLINE __attribute__((always_inline)) inline
#else
#define GLMNETPP_STRONG_INLINE inline
#endif
#endif

/*
 * GLMNETPP_NO_INLINE forces no-inline.
 * using __declspec(noinline) on MSVC, noinline on GCC/clang, and otherwise don't specify an attribute.
 */
#ifndef GLMNETPP_NOINLINE
#if defined(_MSC_VER)
#define GLMNETPP_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define GLMNETPP_NOINLINE __attribute__((noinline)) 
#else
#define GLMNETPP_NOINLINE 
#endif
#endif

/*
 * Generates a CRTP interface.
 */
#ifndef GLMNETPP_GENERATE_CRTP
#define GLMNETPP_GENERATE_CRTP(derived_t) \
    derived_t& self() { return static_cast<derived_t&>(*this); } \
    const derived_t& self() const { return static_cast<const derived_t&>(*this); }
#endif

/*
 * Printing utility.
 */
#ifndef PRINT
#define PRINT(t)                                                         \
 (std::cerr << std::setprecision(18) << __LINE__ << ": " << #t << '\n' \
            << t << "\n"                                              \
            << std::endl)
#endif
