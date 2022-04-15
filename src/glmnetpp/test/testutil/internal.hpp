#pragma once

namespace glmnetpp {

struct InternalParams
{   
    static constexpr double sml = 1e-5;
    static constexpr double eps = 1e-6;
    static constexpr double big = 9.9e35;
    static constexpr int mnlam = 5;
    static constexpr double rsqmax = 0.999;
    static constexpr double pmin = 1e-9;
    static constexpr double exmx = 250.0;
    static constexpr int itrace = 0;
    static constexpr double bnorm_thr = 1e-10;
    static constexpr int bnorm_mxit = 100;
};

} // namespace glmnetpp
