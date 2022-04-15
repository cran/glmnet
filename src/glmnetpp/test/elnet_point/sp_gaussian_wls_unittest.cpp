#include <testutil/base_fixture.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_wls.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_wls.hpp>
#include <testutil/mapped_sparse_matrix_wrapper.hpp>
#include <testutil/translation/wls.hpp>
#include <testutil/data_util.hpp>

namespace glmnetpp {

struct SpWLSPack 
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

    SpWLSPack(double _aint,
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
        using internal_t = SpElnetPointInternal<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::wls,
            double, int, int>;
        using elnet_point_t = SpElnetPoint<util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::wls,
            internal_t>;
        elnet_point_t elnet_point(
                alm0, almc, alpha, X, r, xm, xs, xv, w, intr, ju, vp, 
                cl, nx, thr, maxit, a, aint, g,
                ia, iy, iz, mm, nino, rsqc, nlp);
        elnet_point.fit(m, jerr);
    }

    void fit_old() {
        Eigen::Map<Eigen::VectorXd> xm_wrap(xm.data(), xm.size());
        Eigen::Map<Eigen::VectorXd> xs_wrap(xs.data(), xs.size());
        auto X_wrap = make_mapped_sparse_matrix_wrapper(X, xm_wrap, xs_wrap);
        transl::wls(alm0, almc, alpha, m, X_wrap, r, xv, w, intr, ju, vp, 
            cl, nx, thr, maxit, a, aint, g,
            ia, iy, iz, mm, nino, rsqc, nlp, jerr);
    }
};

struct sp_wls_fixture : 
    base_fixture,
    ::testing::WithParamInterface<
        std::tuple<int, int, int, int, int>
    >
{
protected:

    void test_wls_pack(const SpWLSPack& actual,
                       const SpWLSPack& expected)
    {
        expect_near_vec(actual.r, expected.r, 1e-15);
        expect_near_vec(actual.a, expected.a, 1e-15);
        EXPECT_NEAR(actual.aint, expected.aint, 1e-15);
        expect_near_vec(actual.g, expected.g, 1e-15);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.iy, expected.iy);
        EXPECT_EQ(actual.iz, expected.iz);
        expect_eq_vec(actual.mm, expected.mm);
        EXPECT_EQ(actual.nino, expected.nino);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_NEAR(actual.rsqc, expected.rsqc, actual.rsqc * 1e-15);
    }
        
};

TEST_P(sp_wls_fixture, sp_wls_test)
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
    Eigen::VectorXd xv(p); xv.setOnes();
    Eigen::VectorXi ia(p); ia.setZero();
    Eigen::VectorXi iy(p); iy.setZero();
    Eigen::VectorXi mm(p); mm.setZero();
    Eigen::VectorXd xm(p);
    Eigen::VectorXd xs(p); 
    for (int i = 0; i < xs.size(); ++i) {
        xm(i) = X.col(i).sum() / n;
        xs(i) = X.col(i).cwiseProduct(X.col(i)).sum() / n - xm(i)*xm(i);
    }

    for (int i = 0; i < std::min(p, 10); ++i) {
        ia(i) = mm(i) = i+1;
    }

    SpWLSPack actual(
            aint, nx, X, cl, vp, weights, inclusion, beta, y,
            a, g, r, xm, xs, xv, ia, iy, mm,
            alm0, almc, alpha, intr, maxit);
    SpWLSPack expected(actual);
    run(actual, expected, 0, 1);
    test_wls_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpWLSSuite, sp_wls_fixture,

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
