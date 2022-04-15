#include <legacy/legacy.h>
#include <testutil/translation/multlognetn.hpp>
#include <translation/lognet_base_fixture.hpp>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct MultLognetnPack: LognetBasePack
{
    const Eigen::MatrixXd& y;
    const Eigen::MatrixXd& X;
    const Eigen::VectorXd& xv;
    const int nc;
    Eigen::MatrixXd g;  // g is offset
    Eigen::MatrixXd a0;
    Eigen::VectorXd a;  // flattened 3-d array

    MultLognetnPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            int _intr,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : LognetBasePack(_maxit, _nx, _ne, _nlam, 2, true /*dummy*/, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , xv(_xv)
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
        transl::multlognetn<float, false>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, intr, maxit, xv, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        int nc = y.cols();
        multlognetn_(const_cast<double*>(&alpha), &no, &ni, &nc,
                const_cast<double*>(X.data()), 
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
                &lmu,
                a0.data(), a.data(), m.data(), kin.data(), &dev0, dev.data(), alm.data(),
                &nlp, &jerr);
    }
};

struct multlognetn_fixture
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
        X = dgen.make_X(n, p);
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

        Eigen::VectorXd xm(p), xs(p);
        xv.setZero(p);
        MultLStandardize1::eval(X, w, ju, isd, intr, xm, xs, xv);
    }

protected:
    using base_t = lognet_base_fixture;
    Eigen::MatrixXd X, cl, y;
    Eigen::VectorXd xv, ulam, vp, w;
    Eigen::VectorXi ju;
    int nx, ne, maxit, nlam, intr;
    double alpha, flmin;

    void check_pack(const MultLognetnPack& actual,
                    const MultLognetnPack& expected)
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

TEST_P(multlognetn_fixture, multlognetn_test)
{
    MultLognetnPack actual(
            maxit, nx, ne, nlam, alpha, flmin, intr,
            X, y, w, ulam, xv, vp, cl, ju);
    MultLognetnPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultLognetnSuite, multlognetn_fixture,
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
