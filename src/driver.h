#ifndef DRIVER_H
#define DRIVER_H
#include <glmnetpp_bits/util/exceptions.hpp>
#include <new>

template <class F>
inline void run(F f, int& jerr)
{
    try {
        f();
    }
    catch (const std::bad_alloc&) {
        jerr = glmnetpp::util::bad_alloc_error().err_code(); 
    }
    catch (const std::exception&) {
        jerr = 10001;
    }
}

#endif
