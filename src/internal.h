#pragma once

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
    static double bnorm_thr;
    static int bnorm_mxit;
    static double epsnr;
    static int mxitnr;

    // Formerly stored in the R-level .glmnet_internal env. Centralised
    // here so glmnet.control() reads/writes all 17 algorithm-control
    // parameters through one C++ aggregator. NA_INTEGER on dfmax/pmax
    // means "not overridden -- caller resolves data-dependent default".
    static double thresh;
    static int    maxit;
    static int    dfmax;
    static int    pmax;
};
