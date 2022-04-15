#include <legacy/legacy.h>
#include <testutil/translation/multsplognetn.hpp>
#include <translation/lognet_base_fixture.hpp>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct MultSpLognetnPack: LognetBasePack
{
    const Eigen::MatrixXd& y;
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& xv;
    const Eigen::VectorXd& xm;
    const Eigen::VectorXd& xs;
    const int nc;
    Eigen::MatrixXd g;  // g is offset
    Eigen::MatrixXd a0;
    Eigen::VectorXd a;  // flattened 3-d array

    MultSpLognetnPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            int _intr,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _xm,
            const Eigen::VectorXd& _xs,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : LognetBasePack(_maxit, _nx, _ne, _nlam, 2, true /*dummy*/, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , xv(_xv)
        , xm(_xm)
        , xs(_xs)
        , nc(_y.cols())
        , g(_y.rows(), _y.cols())
        , a0(_y.cols(), _nlam)
        , a(_nx * _y.cols() * _nlam)
    {
        g.setZero(); 
        a0.setZero();
        a.setZero();
    }

    void fit() override
    {
        // don't apply suggested changes
        transl::multsplognetn<float, false>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, intr, maxit, xv, xm, xs, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        int nc = y.cols();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        multsprlognetn_(const_cast<double*>(&alpha), &no, &ni, &nc,
                const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(), 
                const_cast<double*>(y.data()),
                g.data(),
                const_cast<double*>(w.data()),
                const_cast<int*>(ju.data()), 
                const_cast<double*>(vp.data()), 
                const_cast<double*>(cl.data()),
                const_cast<int*>(&ne), 
                const_cast<int*>(&nx), 
                const_cast<int*>(&nlam),
                const_cast<double*>(&flmin), 
                const_cast<double*>(ulam.data()), 
                const_cast<double*>(&thr), 
                const_cast<int*>(&intr),
                const_cast<int*>(&maxit), 
                const_cast<double*>(xv.data()),
                const_cast<double*>(xm.data()),
                const_cast<double*>(xs.data()),
                &lmu,
                a0.data(), a.data(), m.data(), kin.data(), &dev0, dev.data(), alm.data(),
                &nlp, &jerr);
    }
};

struct multsplognetn_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        int isd;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X_sparse(n, p);
        auto beta = dgen.make_beta(p);
        y.resize(n, 2);
        y.col(0) = dgen.make_y(X, beta);
        y.col(0).array() =  (y.col(0).array() > y.col(0).mean()).template cast<double>();
        y.col(1).array() = 1.-y.col(0).array();
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);

        xm.setZero(p);
        xs.setZero(p);
        xv.setZero(p);
        MultSpLStandardize2::eval(X, w, ju, isd, intr, xm, xs, xv);
    }

protected:
    using base_t = lognet_base_fixture;
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd cl, y;
    Eigen::VectorXd xv, ulam, vp, w, xm, xs;
    Eigen::VectorXi ju;
    int nx, ne, maxit, nlam, intr;
    double alpha, flmin;

    void check_pack(const MultSpLognetnPack& actual,
                    const MultSpLognetnPack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_FLOAT_EQ(actual.dev0, expected.dev0);

        expect_eq_vec(actual.kin, expected.kin);
        expect_eq_vec(actual.m, expected.m);

        expect_near_vec(actual.dev, expected.dev, 1e-7);
        expect_float_eq_vec(actual.alm, expected.alm);
        expect_near_mat(actual.g, expected.g, 1e-7);

        expect_near_mat(actual.a0, expected.a0, 1e-7);
        expect_near_vec(actual.a, expected.a, 1e-7);
    }
};

TEST_P(multsplognetn_fixture, multsplognetn_test)
{
    MultSpLognetnPack actual(
            maxit, nx, ne, nlam, alpha, flmin, intr,
            X, y, w, ulam, xv, xm, xs, vp, cl, ju);
    MultSpLognetnPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultSpLognetnSuite, multsplognetn_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5),       // flmin
        testing::Bool(),                      // isd
        testing::Bool()                       // intr
        )
);

} // namespace glmnetpp
