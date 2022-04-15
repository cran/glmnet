#include <testutil/base_fixture.hpp>
#include <testutil/translation/wls.hpp>
#include <testutil/data_util.hpp>
#include <testutil/mapped_sparse_matrix_wrapper.hpp>
#include <testutil/thread.hpp>

namespace glmnetpp {

struct WLSDSCheckPack 
{
    // these values either really won't impact the algorithm
    // or just needs to be default-initialized
    int jerr = 0;
    int nlp = 0;
    int nino = 0;
    int m = 1242;
    double rsqc = 0;
    double thr = 1e-14;
    int iz = 1;

    // initialized upon construction
    double aint;
    int nx;
    const Eigen::SparseMatrix<double>& X;
    const Eigen::MatrixXd& cl;
    const Eigen::VectorXd& vp, w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd beta, y, a, g, r, xm, xs, xv;
    Eigen::VectorXi ia ,iy, mm;
    double alm0;
    double almc;
    double alpha;
    int intr;
    int maxit;

    WLSDSCheckPack(double _aint,
            int _nx,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXd& _vp,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            const Eigen::VectorXd& _beta,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _a,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _r,
            const Eigen::VectorXd& _xm,
            const Eigen::VectorXd& _xs,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXi& _ia,
            const Eigen::VectorXi& _iy,
            const Eigen::VectorXi& _mm,
            double _alm0,
            double _almc,
            double _alpha,
            int _intr,
            int _maxit)
        : aint(_aint), nx(_nx), X(_X), cl(_cl), vp(_vp), w(_w), ju(_ju)
        , beta(_beta), y(_y), a(_a), g(_g), r(_r), xm(_xm), xs(_xs), xv(_xv)
        , ia(_ia), iy(_iy), mm(_mm)
        , alm0(_alm0)
        , almc(_almc)
        , alpha(_alpha)
        , intr(_intr)
        , maxit(_maxit)
    {}

    void fit() {
        Eigen::MatrixXd X_d = X;
        for (int i = 0; i < X_d.cols(); ++i) {
            X_d.col(i).array() -= xm(i);
            X_d.col(i) /= xs(i);
        }
        transl::wls(alm0, almc, alpha, m, X_d, r, xv, w, intr, ju, vp, 
            cl, nx, thr, maxit, a, aint, g,
            ia, iy, iz, mm, nino, rsqc, nlp, jerr);
    }

    void fit_old() {
        Eigen::Map<const Eigen::SparseMatrix<double> > X_m(
                X.rows(), X.cols(), X.nonZeros(),
                X.outerIndexPtr(), X.innerIndexPtr(),
                X.valuePtr(), X.innerNonZeroPtr());
        Eigen::Map<Eigen::VectorXd> xm_m(xm.data(), xm.rows(), xm.cols());
        Eigen::Map<Eigen::VectorXd> xs_m(xs.data(), xs.rows(), xs.cols());
        auto X_w = make_mapped_sparse_matrix_wrapper(X_m, xm_m, xs_m);
        transl::wls(alm0, almc, alpha, m, X_w, r, xv, w, intr, ju, vp, 
            cl, nx, thr, maxit, a, aint, g,
            ia, iy, iz, mm, nino, rsqc, nlp, jerr);
    }
};

struct wls_dscheck_fixture : 
    base_fixture,
    ::testing::WithParamInterface<
        std::tuple<int, int, int, int, int>
    >
{
protected:
    void test_wls_dscheck_pack(const WLSDSCheckPack& actual,
                       const WLSDSCheckPack& expected)
    {
        expect_near_vec(actual.r, expected.r, 1e-14);
        expect_near_vec(actual.a, expected.a, 1e-14);
        EXPECT_NEAR(actual.aint, expected.aint, 1e-14);
        expect_near_vec(actual.g, expected.g, actual.g.array().maxCoeff() * 2e-15);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.iy, expected.iy);
        EXPECT_EQ(actual.iz, expected.iz);
        expect_eq_vec(actual.mm, expected.mm);
        EXPECT_EQ(actual.nino, expected.nino);
        EXPECT_EQ(actual.nlp, expected.nlp);

        // only check error code if it wasn't because of max-active-reached.
        // WLS doesn't log the error code in this case...
        if (actual.jerr > util::max_active_reached_error().err_code(0) ||
            actual.jerr < util::max_active_reached_error().err_code(actual.m)) {
            EXPECT_EQ(actual.jerr, expected.jerr);
        }

        EXPECT_NEAR(actual.rsqc, expected.rsqc, actual.rsqc * 1e-14);
    }
        
};

TEST_P(wls_dscheck_fixture, wls_dscheck_test)
{
    double alm0 = 0.002;
    double almc = 0.001;
    double alpha = 0.9;

    int seed, n, p;
    bool intr, maxit;
    std::tie(seed, n, p, intr, maxit) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X_sparse(n, p);
    auto cl = dgen.make_cl(p);
    auto vp = dgen.make_vp(p);
    auto weights = dgen.make_w(n);
    auto inclusion = dgen.make_ju(p);
    auto beta = dgen.make_beta(p);
    auto y = dgen.make_y(X, beta);
    double aint = 0;
    auto nx = dgen.make_nx(p);

    Eigen::VectorXd a(p); a.setZero();
    Eigen::VectorXd g = (X.transpose() * y).array().abs().matrix();
    Eigen::VectorXd r = y;

    Eigen::MatrixXd X_d = X;
    Eigen::VectorXd xm(p);
    Eigen::VectorXd xs(p); 
    for (int i = 0; i < xs.size(); ++i) {
        if (intr) {
            xm(i) = X_d.col(i).cwiseProduct(weights).mean();
        } else {
            xm(i) = 0.0;
        }
        xs(i) = X_d.col(i).cwiseProduct(X_d.col(i)).cwiseProduct(weights).mean();
        xs(i) -= xm(i) * xm(i);
        xs(i) = std::sqrt(xs(i));
    }
    Eigen::VectorXd xv(p); xv.setOnes();
    Eigen::VectorXi ia(p); ia.setZero();
    Eigen::VectorXi iy(p); iy.setZero();
    Eigen::VectorXi mm(p); mm.setZero();

    for (int i = 0; i < std::min(p, 10); ++i) {
        ia(i) = mm(i) = i+1;
    }

    WLSDSCheckPack actual(
            aint, nx, X, cl, vp, weights, inclusion, beta, y,
            a, g, r, xm, xs, xv, ia, iy, mm,
            alm0, almc, alpha, intr, maxit);
    WLSDSCheckPack expected(actual);
    run(actual, expected, 0, 1);
    test_wls_dscheck_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    WLSDSCheckSuite, wls_dscheck_fixture,

    // combination of inputs: (seed, n, p, intr, maxit)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(5, 10, 20, 30, 40),
        testing::Values(0, 1),
        testing::Values(1, 50, 100)
        )
);

} // namespace glmnetpp
