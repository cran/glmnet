#pragma once
#include <exception>
#include <stdexcept>

namespace glmnetpp {
namespace util {

struct elnet_error : std::exception
{
    virtual int err_code(int) const =0;
};

struct maxit_reached_error : elnet_error
{
    int err_code(int m) const override { return -m-1; } 
};

struct max_active_reached_error : elnet_error 
{
    int err_code(int m) const override { return -10000-m-1; }
};

struct bad_alloc_error : elnet_error
{
    int err_code(int=0) const override { return 1; }
};

} // namespace util
} // namespace glmnetpp
