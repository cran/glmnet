#include <glmnetpp>
#include "tools/bazel_data_util.hpp"
#include <testutil/translation/lognet.hpp>

int main(int argc, char** argv)
{
    using namespace glmnetpp;

    std::string data_path = "bug_reports/data/2022_11_07/";

    Eigen::MatrixXd X = load_csv(argv[0], data_path + "X.csv");
    Eigen::MatrixXd y = load_csv(argv[0], data_path + "y.csv");
    size_t n = X.rows();
    size_t p = X.cols();
    size_t nc = y.cols();

    double alpha = 1.0;
    Eigen::VectorXi jd(1);
    jd.setZero();
    Eigen::MatrixXd offset(n, nc);
    offset.setZero();
    Eigen::VectorXd vp(p);
    vp.setOnes();
    Eigen::MatrixXd cl(2, p);
    cl.row(0).array() = -9.9e35;
    cl.row(1).array() = 9.9e35;
    int ne = 129;
    int nx = 128;
    int nlam = 100;
    double flmin = 1e-4;
    Eigen::VectorXd ulam(1); 
    ulam[0] = 0;
    double thr = 1e-7;
    bool isd = true;
    bool intr = true;
    int maxit = 1e5;
    int kopt = 0;

    int lmu;
    Eigen::MatrixXd a0(nc, nlam);
    a0.setZero();
    Eigen::VectorXd ca(nx * nlam * nc);
    Eigen::VectorXi ia(nx);
    Eigen::VectorXi nin(nlam);
    double dev0;
    Eigen::VectorXd dev(nlam);
    Eigen::VectorXd alm(nlam);
    int nlp;
    int jerr;
    
    transl::lognet<double, true>(
        alpha, X, y, offset, jd, vp, cl, ne, nx, nlam, flmin,
        ulam, thr, isd, intr, maxit, kopt, lmu, a0, ca, ia,
        nin, dev0, dev, alm, nlp, jerr
    );
    
    return 0;
}
