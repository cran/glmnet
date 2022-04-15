#include "internal.h"
#include <cstddef>
#include <Rcpp.h>
#include <RcppEigen.h>
// ORDER MATTERS
#include <R.h>
#include <Rinternals.h>

extern "C" {

void F77_SUB(chg_fract_dev)(double*); 
void F77_SUB(chg_dev_max)(double*); 
void F77_SUB(chg_min_flmin)(double*); 
void F77_SUB(chg_big)(double*); 
void F77_SUB(chg_min_lambdas)(int*); 
void F77_SUB(chg_min_null_prob)(double*); 
void F77_SUB(chg_max_exp)(double*); 
void F77_SUB(chg_itrace)(int*); 
void F77_SUB(chg_bnorm)(double*, int*); 
void F77_SUB(chg_epsnr)(double*); 
void F77_SUB(chg_mxitnr)(int*); 

} // end extern "C"

double InternalParams::sml = 1e-5;
double InternalParams::eps = 1e-6;
double InternalParams::big = 9.9e35;
int InternalParams::mnlam = 5;
double InternalParams::rsqmax = 0.999;
double InternalParams::pmin = 1e-9;
double InternalParams::exmx = 250.0;
int InternalParams::itrace = 0;
double InternalParams::bnorm_thr = 1e-10;
int InternalParams::bnorm_mxit = 100;
double InternalParams::epsnr = 1e-6;
int InternalParams::mxitnr = 25;

using namespace Rcpp;

// [[Rcpp::export]]
List get_int_parms(double& fdev,
                   double& eps,
                   double& big,
                   int& mnlam,
                   double& devmax,
                   double& pmin,
                   double& exmx,
                   int& itrace)
{
    fdev = InternalParams::sml; 
    eps = InternalParams::eps; 
    big = InternalParams::big; 
    mnlam = InternalParams::mnlam; 
    devmax = InternalParams::rsqmax;
    pmin = InternalParams::pmin; 
    exmx = InternalParams::exmx; 
    itrace = InternalParams::itrace;
    return List::create(
            Named("fdev")=fdev,
            Named("eps")=eps,
            Named("big")=big,
            Named("mnlam")=mnlam,
            Named("devmax")=devmax,
            Named("pmin")=pmin,
            Named("exmx")=exmx,
            Named("itrace")=itrace);
}

// [[Rcpp::export]]
List get_int_parms2(double& epsnr, int& mxitnr)
{
    epsnr = InternalParams::epsnr;
    mxitnr = InternalParams::mxitnr;
    return List::create(
            Named("epsnr")=epsnr,
            Named("mxitnr")=mxitnr);
}

// [[Rcpp::export]]
void chg_fract_dev(double arg) { /*TODO*/ F77_SUB(chg_fract_dev)(&arg); InternalParams::sml = arg; }

// [[Rcpp::export]]
void chg_dev_max(double arg) { /*TODO*/ F77_SUB(chg_dev_max)(&arg); InternalParams::rsqmax = arg; }

// [[Rcpp::export]]
void chg_min_flmin(double arg) { /*TODO*/ F77_SUB(chg_min_flmin)(&arg); InternalParams::eps = arg; }

// [[Rcpp::export]]
void chg_big(double arg) { /*TODO*/ F77_SUB(chg_big)(&arg); InternalParams::big = arg; }

// [[Rcpp::export]]
void chg_min_lambdas(int irg) { /*TODO*/ F77_SUB(chg_min_lambdas)(&irg); InternalParams::mnlam = irg; }

// [[Rcpp::export]]
void chg_min_null_prob(double arg) { /*TODO*/ F77_SUB(chg_min_null_prob)(&arg); InternalParams::pmin = arg; }

// [[Rcpp::export]]
void chg_max_exp(double arg) { /*TODO*/ F77_SUB(chg_max_exp)(&arg); InternalParams::exmx = arg; }

// [[Rcpp::export]]
void chg_itrace(int irg) { /*TODO*/ F77_SUB(chg_itrace)(&irg); InternalParams::itrace = irg; }

// [[Rcpp::export]]
void chg_bnorm(double arg, int irg) { 
    /*TODO*/
    F77_SUB(chg_bnorm)(&arg, &irg);
    InternalParams::bnorm_thr = arg; 
    InternalParams::bnorm_mxit = irg; 
}

// [[Rcpp::export]]
List get_bnorm(double& prec, int& mxit) {
    prec = InternalParams::bnorm_thr; 
    mxit = InternalParams::bnorm_mxit;
    return List::create(
            Named("prec")=prec,
            Named("mxit")=mxit);
}

// [[Rcpp::export]]
void chg_epsnr(double arg) { /*TODO*/ F77_SUB(chg_epsnr)(&arg); InternalParams::epsnr = arg; }

// [[Rcpp::export]]
void chg_mxitnr(int irg) { /*TODO*/ F77_SUB(chg_mxitnr)(&irg); InternalParams::mxitnr = irg; }

