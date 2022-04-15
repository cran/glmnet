#include <testutil/base_fixture.hpp>
#include <testutil/translation/wls.hpp>
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <testutil/thread.hpp>

namespace glmnetpp {

struct WLSPack 
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
    const Eigen::MatrixXd& X, cl;
    const Eigen::VectorXd& vp, w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd beta, y, a, g, r, xv;
    Eigen::VectorXi ia ,iy, mm;
    double alm0;
    double almc;
    double alpha;
    int intr;
    int maxit;

    WLSPack(double _aint,
            int _nx,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXd& _vp,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            const Eigen::VectorXd& _beta,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _a,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _r,
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
        , beta(_beta), y(_y), a(_a), g(_g), r(_r), xv(_xv)
        , ia(_ia), iy(_iy), mm(_mm)
        , alm0(_alm0)
        , almc(_almc)
        , alpha(_alpha)
        , intr(_intr)
        , maxit(_maxit)
    {}

    void fit() {
        transl::wls(alm0, almc, alpha, m, X, r, xv, w, intr, ju, vp, 
            cl, nx, thr, maxit, a, aint, g,
            ia, iy, iz, mm, nino, rsqc, nlp, jerr);
    }

    void fit_old() {
        int n = X.rows();
        int p = X.cols();
        wls_(&alm0, &almc, &alpha, &m, &n, &p, 
             const_cast<double*>(X.data()), r.data(), 
             const_cast<double*>(w.data()), &intr, 
             const_cast<int*>(ju.data()), 
             const_cast<double*>(vp.data()), 
             const_cast<double*>(cl.data()), 
             &nx, &thr, &maxit, a.data(), &aint, g.data(),
             ia.data(), iy.data(), &iz, mm.data(), &nino, &rsqc, &nlp, &jerr);
    }
};

struct wls_fixture : 
    base_fixture,
    ::testing::WithParamInterface<
        std::tuple<int, int, int, int, int>
    >
{
protected:
    void test_wls_pack(const WLSPack& actual,
                       const WLSPack& expected)
    {
        expect_near_vec(actual.r, expected.r, 1e-14);
        expect_near_vec(actual.a, expected.a, 1e-14);
        EXPECT_NEAR(actual.aint, expected.aint, 1e-14);
        expect_near_vec(actual.g, expected.g, 1e-14);
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

TEST_P(wls_fixture, wls_test)
{
    double alm0 = 0.002;
    double almc = 0.001;
    double alpha = 0.9;

    int seed, n, p;
    bool intr, maxit;
    std::tie(seed, n, p, intr, maxit) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X(n, p);
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
    Eigen::VectorXd xv(p); xv.setOnes();
    Eigen::VectorXi ia(p); ia.setZero();
    Eigen::VectorXi iy(p); iy.setZero();
    Eigen::VectorXi mm(p); mm.setZero();

    for (int i = 0; i < std::min(p, 10); ++i) {
        ia(i) = mm(i) = i+1;
    }

    WLSPack actual(
            aint, nx, X, cl, vp, weights, inclusion, beta, y,
            a, g, r, xv, ia, iy, mm,
            alm0, almc, alpha, intr, maxit);
    WLSPack expected(actual);
    run(actual, expected, 0, 1);
    test_wls_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    WLSSuite, wls_fixture,

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
