#pragma once

namespace glmnetpp {

struct InternalParams
{
    static double sml;
    static double eps;
    static double big;
    static int mnlam;
    static double rsqmax;
    static double pmin;
    static double exmx;
    static int itrace;
};

void get_int_parms(double& sml,
                   double& eps,
                   double& big,
                   int& mnlam,
                   double& rsqmax,
                   double& pmin,
                   double& exmx,
                   int& itrace);

void chg_fract_dev(double arg);
void chg_min_flmin(double arg); 
void chg_dev_max(double arg);
void chg_big(double arg);
void chg_min_lambdas(int irg);
void chg_min_null_prob(double arg);
void chg_max_exp(double arg);
void chg_itrace(int irg);

} // namespace glmnetpp
