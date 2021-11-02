#include <glmnetpp_bits/internal.hpp>

namespace glmnetpp {

double InternalParams::sml = 1e-5;
double InternalParams::eps = 1e-6;
double InternalParams::big = 9.9e35;
int InternalParams::mnlam = 5;
double InternalParams::rsqmax = 0.999;
double InternalParams::pmin = 1e-9;
double InternalParams::exmx = 250.0;
int InternalParams::itrace = 0;

// TODO: this interface is kinda terrible,
// but need it to be compatible with current R interface.

void get_int_parms(double& sml,
                   double& eps,
                   double& big,
                   int& mnlam,
                   double& rsqmax,
                   double& pmin,
                   double& exmx,
                   int& itrace)
{
    sml = InternalParams::sml; 
    eps = InternalParams::eps; 
    big = InternalParams::big; 
    mnlam = InternalParams::mnlam; 
    rsqmax = InternalParams::rsqmax;
    pmin = InternalParams::pmin; 
    exmx = InternalParams::exmx; 
    itrace = InternalParams::itrace;
}

void chg_fract_dev(double arg) { InternalParams::sml = arg; }
void chg_min_flmin(double arg) { InternalParams::eps = arg; }
void chg_dev_max(double arg) { InternalParams::rsqmax = arg; }
void chg_big(double arg) { InternalParams::big = arg; }
void chg_min_lambdas(int irg) { InternalParams::mnlam = irg; }
void chg_min_null_prob(double arg) { InternalParams::pmin = arg; }
void chg_max_exp(double arg) { InternalParams::exmx = arg; }
void chg_itrace(int irg) { InternalParams::itrace = irg; }

} // namespace glmnetpp
