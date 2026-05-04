// Aggregated getter and setter for glmnet's algorithm-control parameters.
// One Rcpp::List in/out, one .Call entry each. The R-side glmnet.control()
// (see R/glmnet.control.R) is a thin wrapper around these two.
//
// Field map (canonical order matches glmnet.control()'s historical return):
//
//   R name      C++ storage                   type
//   -------     ----------------              ----
//   fdev        InternalParams::sml           double
//   eps         InternalParams::eps           double
//   big         InternalParams::big           double
//   mnlam       InternalParams::mnlam         int
//   devmax      InternalParams::rsqmax        double
//   pmin        InternalParams::pmin          double
//   exmx        InternalParams::exmx          double
//   itrace      InternalParams::itrace        int
//   prec        InternalParams::bnorm_thr     double
//   mxit        InternalParams::bnorm_mxit    int
//   epsnr       InternalParams::epsnr         double
//   mxitnr      InternalParams::mxitnr        int
//   thresh      InternalParams::thresh        double
//   maxit       InternalParams::maxit         int
//   dfmax       InternalParams::dfmax         int (NA_INTEGER == NULL)
//   pmax        InternalParams::pmax          int (NA_INTEGER == NULL)
//   trace.it    InternalParams::itrace        int (alias of itrace)
//
// Unknown keys passed to glmnet_control_set are silently ignored (per
// design): no R-side validation, no error.

#include "internal.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List glmnet_control_get() {
    return List::create(
        Named("fdev")     = InternalParams::sml,
        Named("eps")      = InternalParams::eps,
        Named("big")      = InternalParams::big,
        Named("mnlam")    = InternalParams::mnlam,
        Named("devmax")   = InternalParams::rsqmax,
        Named("pmin")     = InternalParams::pmin,
        Named("exmx")     = InternalParams::exmx,
        Named("itrace")   = InternalParams::itrace,
        Named("prec")     = InternalParams::bnorm_thr,
        Named("mxit")     = InternalParams::bnorm_mxit,
        Named("epsnr")    = InternalParams::epsnr,
        Named("mxitnr")   = InternalParams::mxitnr,
        Named("thresh")   = InternalParams::thresh,
        Named("maxit")    = InternalParams::maxit,
        Named("dfmax")    = InternalParams::dfmax,
        Named("pmax")     = InternalParams::pmax,
        Named("trace.it") = InternalParams::itrace
    );
}

// [[Rcpp::export]]
void glmnet_control_set(List updates) {
    if (updates.size() == 0) return;
    CharacterVector names = updates.names();
    int n = names.size();
    for (int i = 0; i < n; ++i) {
        std::string key = as<std::string>(names[i]);
        SEXP v = updates[i];

        // dfmax/pmax accept NULL (means "unset, resolve at call site").
        if (key == "dfmax") {
            InternalParams::dfmax = Rf_isNull(v) ? NA_INTEGER : as<int>(v);
            continue;
        }
        if (key == "pmax") {
            InternalParams::pmax  = Rf_isNull(v) ? NA_INTEGER : as<int>(v);
            continue;
        }

        // For all other fields, NULL is silently ignored (no-op).
        if (Rf_isNull(v)) continue;

        if      (key == "fdev")     InternalParams::sml        = as<double>(v);
        else if (key == "eps")      InternalParams::eps        = as<double>(v);
        else if (key == "big")      InternalParams::big        = as<double>(v);
        else if (key == "mnlam")    InternalParams::mnlam      = as<int>(v);
        else if (key == "devmax")   InternalParams::rsqmax     = as<double>(v);
        else if (key == "pmin")     InternalParams::pmin       = as<double>(v);
        else if (key == "exmx")     InternalParams::exmx       = as<double>(v);
        else if (key == "itrace")   InternalParams::itrace     = as<int>(v);
        else if (key == "prec")     InternalParams::bnorm_thr  = as<double>(v);
        else if (key == "mxit")     InternalParams::bnorm_mxit = as<int>(v);
        else if (key == "epsnr")    InternalParams::epsnr      = as<double>(v);
        else if (key == "mxitnr")   InternalParams::mxitnr     = as<int>(v);
        else if (key == "thresh")   InternalParams::thresh     = as<double>(v);
        else if (key == "maxit")    InternalParams::maxit      = as<int>(v);
        else if (key == "trace.it") InternalParams::itrace     = as<int>(v);
        // unknown keys silently ignored
    }
}
