#include <cstddef>
#include <RcppEigen.h>
#include <glmnetpp>
#include <R.h>
#include <Rinternals.h>
#include "internal_params.h"

using namespace Rcpp;
using namespace glmnetpp;

void setpb_cpp(SEXP, int);

extern "C" 
{
void F77_SUB(get_int_parms)(
			    double *sml,
			    double *eps,
			    double *big,
			    int *mnlam,
			    double *rsqmax,
			    double *pmin,
			    double *exmx,
			    int *itrace
			    );
}

// Elnet for dense X.
// [[Rcpp::export]]
List elnet_exp(
    int ka,
    double parm,
    Eigen::MatrixXd x,
    Eigen::VectorXd y,
    Eigen::VectorXd w,          // TODO: figure out if we should allow updating (safe choice is to copy)
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,
    int ne,
    int nx,
    int nlam,
    double flmin,
    const Eigen::Map<Eigen::VectorXd> ulam,
    double thr,
    int isd,
    int intr,
    int maxit,
    SEXP pb,
    int lmu,
    Eigen::Map<Eigen::VectorXd> a0,
    Eigen::Map<Eigen::MatrixXd> ca,
    Eigen::Map<Eigen::VectorXi> ia,
    Eigen::Map<Eigen::VectorXi> nin,
    Eigen::Map<Eigen::VectorXd> rsq,
    Eigen::Map<Eigen::VectorXd> alm,
    int nlp,
    int jerr
    )
{
    using elnet_driver_t = ElnetDriver<util::glm_type::gaussian>;
    elnet_driver_t driver;
    InternalParamsExp int_params;
    F77_SUB(get_int_parms)(&int_params.sml,
                           &int_params.eps,
                           &int_params.big,
                           &int_params.mnlam,
                           &int_params.rsqmax,
                           &int_params.pmin,
                           &int_params.exmx,
                           &int_params.itrace);
    try {
        driver.fit(
                ka == 2, parm, x, y, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd == 1, intr == 1, maxit, 
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr, 
                [=](int v) {setpb_cpp(pb, v);}, int_params);
    }
    catch (const std::bad_alloc&) {
        jerr = util::bad_alloc_error().err_code(); 
    }

    return List::create(
            Named("a0")=a0,
            Named("nin")=nin,
            Named("alm")=alm,
            Named("ca")=ca,
            Named("ia")=ia,
            Named("lmu")=lmu,
            Named("rsq")=rsq,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}

// Elnet for sparse X.
// [[Rcpp::export]]
List spelnet_exp(
    int ka,
    double parm,
    const Eigen::Map<Eigen::SparseMatrix<double> > x,
    Eigen::VectorXd y,
    Eigen::VectorXd w,
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,
    int ne,
    int nx,
    int nlam,
    double flmin,
    const Eigen::Map<Eigen::VectorXd> ulam,
    double thr,
    int isd,
    int intr,
    int maxit,
    SEXP pb,
    int lmu,
    Eigen::Map<Eigen::VectorXd> a0,
    Eigen::Map<Eigen::MatrixXd> ca,
    Eigen::Map<Eigen::VectorXi> ia,
    Eigen::Map<Eigen::VectorXi> nin,
    Eigen::Map<Eigen::VectorXd> rsq,
    Eigen::Map<Eigen::VectorXd> alm,
    int nlp,
    int jerr
    )
{
    using elnet_driver_t = ElnetDriver<util::glm_type::gaussian>;
    elnet_driver_t driver;
    InternalParamsExp int_params;
    F77_SUB(get_int_parms)(&int_params.sml,
                           &int_params.eps,
                           &int_params.big,
                           &int_params.mnlam,
                           &int_params.rsqmax,
                           &int_params.pmin,
                           &int_params.exmx,
                           &int_params.itrace);
    try {
        driver.fit(
                ka == 2, parm, x, y, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd == 1, intr == 1, maxit, 
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr, 
                [=](int v) {setpb_cpp(pb, v);}, int_params);
    }
    catch (const std::bad_alloc&) {
        jerr = util::bad_alloc_error().err_code(); 
    }

    return List::create(
            Named("a0")=a0,
            Named("nin")=nin,
            Named("alm")=alm,
            Named("ca")=ca,
            Named("ia")=ia,
            Named("lmu")=lmu,
            Named("rsq")=rsq,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}
