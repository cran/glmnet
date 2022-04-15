#include <cstddef>
#include <RcppEigen.h>
#include <glmnetpp>
#include <R.h>
#include <Rinternals.h>
#include "driver.h"
#include "internal.h"

using namespace Rcpp;
using namespace glmnetpp;

void setpb_cpp(SEXP, int);

// Gaussian for dense X.
// [[Rcpp::export]]
List elnet_exp(
    int ka,
    double parm,
    Eigen::MatrixXd x,          // TODO: map?
    Eigen::VectorXd y,          // TODO: map?
    Eigen::VectorXd w,          // TODO: figure out if we should allow updating (safe choice is to copy)
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,         // TODO: map?
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
    auto f = [&]() {
        driver.fit(
                ka == 2, parm, x, y, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd == 1, intr == 1, maxit, 
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr, 
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
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

// Gaussian for sparse X.
// [[Rcpp::export]]
List spelnet_exp(
    int ka,
    double parm,
    const Eigen::Map<Eigen::SparseMatrix<double> > x,
    Eigen::VectorXd y, // TODO: map?
    Eigen::VectorXd w, // TODO: map?
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl, // TODO: map?
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
    auto f = [&]() {
        driver.fit(
                ka == 2, parm, x, y, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd == 1, intr == 1, maxit, 
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr, 
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
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

// Binomial/Multinomial for dense X.
// [[Rcpp::export]]
List lognet_exp(
    double parm,
    Eigen::MatrixXd x,          // TODO: map?
    Eigen::MatrixXd y,          // TODO: map?
    Eigen::MatrixXd g,          // TODO: map? 
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,         // TODO: map?
    int ne,
    int nx,
    int nlam,
    double flmin,
    const Eigen::Map<Eigen::VectorXd> ulam,
    double thr,
    int isd,
    int intr,
    int maxit,
    int kopt,
    SEXP pb,
    int lmu,
    Eigen::Map<Eigen::MatrixXd> a0,
    Eigen::Map<Eigen::VectorXd> ca,
    Eigen::Map<Eigen::VectorXi> ia,
    Eigen::Map<Eigen::VectorXi> nin,
    double nulldev,
    Eigen::Map<Eigen::VectorXd> dev,
    Eigen::Map<Eigen::VectorXd> alm,
    int nlp,
    int jerr
    )
{
    using elnet_driver_t = ElnetDriver<util::glm_type::binomial>;
    elnet_driver_t driver;
    auto f = [&]() {
        driver.fit(
                parm, x, y, g, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd, intr, maxit, kopt,
                lmu, a0, ca, ia, nin, nulldev, dev, alm, nlp, jerr,
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
    return List::create(
            Named("a0")=a0,
            Named("nin")=nin,
            Named("alm")=alm,
            Named("ca")=ca,
            Named("ia")=ia,
            Named("lmu")=lmu,
            Named("nulldev")=nulldev,
            Named("dev")=dev,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}

// Lognet for sparse X.
// [[Rcpp::export]]
List splognet_exp(
    double parm,
    const Eigen::Map<Eigen::SparseMatrix<double>> x,
    Eigen::MatrixXd y,          // TODO: map?
    Eigen::MatrixXd g,          // TODO: map? 
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,         // TODO: map?
    int ne,
    int nx,
    int nlam,
    double flmin,
    const Eigen::Map<Eigen::VectorXd> ulam,
    double thr,
    int isd,
    int intr,
    int maxit,
    int kopt,
    SEXP pb,
    int lmu,
    Eigen::Map<Eigen::MatrixXd> a0,
    Eigen::Map<Eigen::VectorXd> ca,
    Eigen::Map<Eigen::VectorXi> ia,
    Eigen::Map<Eigen::VectorXi> nin,
    double nulldev,
    Eigen::Map<Eigen::VectorXd> dev,
    Eigen::Map<Eigen::VectorXd> alm,
    int nlp,
    int jerr
    )
{
    using elnet_driver_t = ElnetDriver<util::glm_type::binomial>;
    elnet_driver_t driver;
    auto f = [&]() {
        driver.fit(
                parm, x, y, g, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd, intr, maxit, kopt,
                lmu, a0, ca, ia, nin, nulldev, dev, alm, nlp, jerr,
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
    return List::create(
            Named("a0")=Eigen::MatrixXd(a0),
            Named("nin")=Eigen::VectorXi(nin),
            Named("alm")=Eigen::VectorXd(alm),
            Named("ca")=Eigen::VectorXd(ca),
            Named("ia")=Eigen::VectorXi(ia),
            Named("lmu")=lmu,
            Named("nulldev")=nulldev,
            Named("dev")=Eigen::VectorXd(dev),
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}

// Poisson for dense X.
// [[Rcpp::export]]
List fishnet_exp(
    double parm,
    Eigen::MatrixXd x,          // TODO: map?
    Eigen::VectorXd y,          // TODO: map?
    Eigen::VectorXd g,          // TODO: map? 
    const Eigen::Map<Eigen::VectorXd> w,
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,         // TODO: map?
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
    double nulldev,
    Eigen::Map<Eigen::VectorXd> dev,
    Eigen::Map<Eigen::VectorXd> alm,
    int nlp,
    int jerr
    )
{
    using elnet_driver_t = ElnetDriver<util::glm_type::poisson>;
    elnet_driver_t driver;
    auto f = [&]() {
        driver.fit(
                parm, x, y, g, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd, intr, maxit,
                lmu, a0, ca, ia, nin, nulldev, dev, alm, nlp, jerr,
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
    return List::create(
            Named("a0")=a0,
            Named("nin")=nin,
            Named("alm")=alm,
            Named("ca")=ca,
            Named("ia")=ia,
            Named("lmu")=lmu,
            Named("nulldev")=nulldev,
            Named("dev")=dev,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}

// Poisson for sparse X.
// [[Rcpp::export]]
List spfishnet_exp(
    double parm,
    const Eigen::Map<Eigen::SparseMatrix<double>> x,
    Eigen::VectorXd y,          // TODO: map?
    Eigen::VectorXd g,          // TODO: map? 
    const Eigen::Map<Eigen::VectorXd> w,
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    Eigen::MatrixXd cl,         // TODO: map?
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
    double nulldev,
    Eigen::Map<Eigen::VectorXd> dev,
    Eigen::Map<Eigen::VectorXd> alm,
    int nlp,
    int jerr
    )
{
    using elnet_driver_t = ElnetDriver<util::glm_type::poisson>;
    elnet_driver_t driver;
    auto f = [&]() {
        driver.fit(
                parm, x, y, g, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd, intr, maxit,
                lmu, a0, ca, ia, nin, nulldev, dev, alm, nlp, jerr,
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
    return List::create(
            Named("a0")=a0,
            Named("nin")=nin,
            Named("alm")=alm,
            Named("ca")=ca,
            Named("ia")=ia,
            Named("lmu")=lmu,
            Named("nulldev")=nulldev,
            Named("dev")=dev,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}

// Multi-response Gaussian for dense X.
// [[Rcpp::export]]
List multelnet_exp(
    double parm,
    Eigen::MatrixXd x,          // TODO: map?
    Eigen::MatrixXd y,          // TODO: map?
    Eigen::VectorXd w,          // TODO: map?
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    const Eigen::Map<Eigen::MatrixXd> cl,
    int ne,
    int nx,
    int nlam,
    double flmin,
    const Eigen::Map<Eigen::VectorXd> ulam,
    double thr,
    int isd,
    int jsd,
    int intr,
    int maxit,
    SEXP pb,
    int lmu,
    Eigen::Map<Eigen::MatrixXd> a0,
    Eigen::Map<Eigen::VectorXd> ca,
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
    auto f = [&]() {
        driver.fit(
                parm, x, y, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd, jsd, intr, maxit,
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr,
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
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

// Multi-response Gaussian for sparse X.
// [[Rcpp::export]]
List multspelnet_exp(
    double parm,
    const Eigen::Map<Eigen::SparseMatrix<double> > x,
    Eigen::MatrixXd y,          // TODO: map?
    Eigen::VectorXd w,          // TODO: map?
    const Eigen::Map<Eigen::VectorXi> jd,
    const Eigen::Map<Eigen::VectorXd> vp,
    const Eigen::Map<Eigen::MatrixXd> cl,
    int ne,
    int nx,
    int nlam,
    double flmin,
    const Eigen::Map<Eigen::VectorXd> ulam,
    double thr,
    int isd,
    int jsd,
    int intr,
    int maxit,
    SEXP pb,
    int lmu,
    Eigen::Map<Eigen::MatrixXd> a0,
    Eigen::Map<Eigen::VectorXd> ca,
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
    auto f = [&]() {
        driver.fit(
                parm, x, y, w, jd, vp, cl, ne, nx, nlam, flmin,
                ulam, thr, isd, jsd, intr, maxit,
                lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr,
                [&](int v) {setpb_cpp(pb, v);}, ::InternalParams());
    };
    run(f, jerr);
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
