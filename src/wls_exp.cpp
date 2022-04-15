#include <cstddef>
#include <RcppEigen.h>
#include <glmnetpp>
#include "driver.h"

using namespace Rcpp;
using namespace glmnetpp;

// WLS for dense X.
// [[Rcpp::export]]
List wls_exp(
    double alm0,
    double almc,
    double alpha,
    int m,
    int no,
    int ni,
    const Eigen::Map<Eigen::MatrixXd> x,
    Eigen::Map<Eigen::VectorXd> r,
    Eigen::Map<Eigen::VectorXd> xv,
    const Eigen::Map<Eigen::VectorXd> v,
    int intr,
    const Eigen::Map<Eigen::VectorXi> ju,
    const Eigen::Map<Eigen::VectorXd> vp,
    const Eigen::Map<Eigen::MatrixXd> cl,
    int nx,
    double thr,
    int maxit,
    Eigen::Map<Eigen::VectorXd> a,
    double aint,
    Eigen::Map<Eigen::VectorXd> g,
    Eigen::Map<Eigen::VectorXi> ia,
    Eigen::Map<Eigen::VectorXi> iy,
    int iz,
    Eigen::Map<Eigen::VectorXi> mm,
    int nino,
    double rsqc,
    int nlp,
    int jerr
    )
{
    using internal_t = ElnetPointInternal<
        util::glm_type::gaussian,
        util::mode_type<util::glm_type::gaussian>::wls,
        double, int, int>;
    using elnet_point_t = ElnetPoint<
        util::glm_type::gaussian,
        util::mode_type<util::glm_type::gaussian>::wls,
        internal_t>;
    auto f = [&]() {
        elnet_point_t elnet_point(
                alm0, almc, alpha, x, r, xv, v, intr, ju, vp,
                cl, nx, thr, maxit, a, aint, g, 
                ia, iy, iz, mm, nino, rsqc, nlp);
        elnet_point.fit(m, jerr);
    };
    run(f, jerr);

    return List::create(
            Named("almc")=almc,
            Named("r")=r,
            Named("xv")=xv,
            Named("ju")=ju,
            Named("vp")=vp,
            Named("cl")=cl,
            Named("nx")=nx,
            Named("a")=a,
            Named("aint")=aint,
            Named("g")=g,
            Named("ia")=ia,
            Named("iy")=iy,
            Named("iz")=iz,
            Named("mm")=mm,
            Named("nino")=nino,
            Named("rsqc")=rsqc,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}

// WLS for sparse X.
// [[Rcpp::export]]
List spwls_exp(
    double alm0,
    double almc,
    double alpha,
    int m,
    int no,
    int ni,
    const Eigen::Map<Eigen::SparseMatrix<double>> x,
    const Eigen::Map<Eigen::VectorXd> xm,
    const Eigen::Map<Eigen::VectorXd> xs,
    Eigen::Map<Eigen::VectorXd> r,
    Eigen::Map<Eigen::VectorXd> xv,
    const Eigen::Map<Eigen::VectorXd> v,
    int intr,
    const Eigen::Map<Eigen::VectorXi> ju,
    const Eigen::Map<Eigen::VectorXd> vp,
    const Eigen::Map<Eigen::MatrixXd> cl,
    int nx,
    double thr,
    int maxit,
    Eigen::Map<Eigen::VectorXd> a,
    double aint,
    Eigen::Map<Eigen::VectorXd> g,
    Eigen::Map<Eigen::VectorXi> ia,
    Eigen::Map<Eigen::VectorXi> iy,
    int iz,
    Eigen::Map<Eigen::VectorXi> mm,
    int nino,
    double rsqc,
    int nlp,
    int jerr
    )
{
    using internal_t = SpElnetPointInternal<
        util::glm_type::gaussian,
        util::mode_type<util::glm_type::gaussian>::wls,
        double, int, int>;
    using elnet_point_t = SpElnetPoint<
        util::glm_type::gaussian,
        util::mode_type<util::glm_type::gaussian>::wls,
        internal_t>;
    auto f = [&]() {
        elnet_point_t elnet_point(
                alm0, almc, alpha, x, r, xm, xs, xv, v, intr, ju, vp,
                cl, nx, thr, maxit, a, aint, g, 
                ia, iy, iz, mm, nino, rsqc, nlp);
        elnet_point.fit(m, jerr);
    };
    run(f, jerr);

    return List::create(
            Named("almc")=almc,
            Named("r")=r,
            Named("xv")=xv,
            Named("ju")=ju,
            Named("vp")=vp,
            Named("cl")=cl,
            Named("nx")=nx,
            Named("a")=a,
            Named("aint")=aint,
            Named("g")=g,
            Named("ia")=ia,
            Named("iy")=iy,
            Named("iz")=iz,
            Named("mm")=mm,
            Named("nino")=nino,
            Named("rsqc")=rsqc,
            Named("nlp")=nlp,
            Named("jerr")=jerr);
}
