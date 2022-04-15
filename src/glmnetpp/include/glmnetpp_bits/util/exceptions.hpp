#pragma once
#include <exception>
#include <stdexcept>

namespace glmnetpp {
namespace util {

// TODO: change interface so that we construct with int and err_code has no args.
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
    int err_code(int m) const override { return -10001-m; }
};

struct bad_alloc_error : elnet_error
{
    int err_code(int=0) const override { return 1; }
};

struct prob_min_reached_error : elnet_error
{
    prob_min_reached_error(int m)
        : m_(m)
    {} 
    int err_code(int=0) const override { return 8001+m_; }

private:
    int m_;
};

struct prob_max_reached_error : elnet_error
{
    prob_max_reached_error(int m)
        : m_(m)
    {}
    int err_code(int=0) const override { return 9001+m_; }

private:
    int m_;
};

struct below_min_variance_error : elnet_error
{
    int err_code(int m) const override { return -20001-m; }
};

struct bnorm_maxit_reached_error : elnet_error
{
    int err_code(int m=0) const override { return 90000; }
};

struct non_positive_penalty_error : elnet_error
{
    int err_code(int m=0) const override { return 10000; }
};

struct negative_response_error : elnet_error
{
    int err_code(int m=0) const override { return 8888; }
};

struct all_excluded_error : elnet_error
{
    int err_code(int m=0) const override { return 7777; }
};

struct non_positive_weight_sum_error : elnet_error
{
    int err_code(int m=0) const override { return 9999; }
};

} // namespace util
} // namespace glmnetpp
