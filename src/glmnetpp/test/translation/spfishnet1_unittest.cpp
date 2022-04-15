#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <testutil/translation/spfishnet1.hpp>
#include <testutil/thread.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <Eigen/SparseCore>

namespace glmnetpp {

struct SpFishnet1Pack
{
    const double beta, flmin, thr = 1e-14;
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& y;
    const Eigen::VectorXd& q;
    const Eigen::VectorXi& ju;
    const Eigen::VectorXd& vp;
    const Eigen::MatrixXd& cl;
    const int ne, nx, nlam, maxit;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& xb;
    const Eigen::VectorXd& xs;
    const bool intr; 

    int lmu = 0, nlp = 0, jerr = 0;
    Eigen::VectorXd g;
    Eigen::VectorXd a0;
    Eigen::MatrixXd ca;
    Eigen::VectorXi ia;
    Eigen::VectorXi kin;
    double dev0 = 0;
    Eigen::VectorXd dev;
    Eigen::VectorXd alm;

    SpFishnet1Pack(
        double _beta, 
        double _flmin, 
        const Eigen::SparseMatrix<double>& _X,
        const Eigen::VectorXd& _y,
        const Eigen::VectorXd& _g,
        const Eigen::VectorXd& _q,
        const Eigen::VectorXi& _ju,
        const Eigen::VectorXd& _vp,
        const Eigen::MatrixXd& _cl,
        int _ne, 
        int _nx, 
        int _nlam, 
        int _maxit,
        const Eigen::VectorXd& _ulam,
        const Eigen::VectorXd& _xb,
        const Eigen::VectorXd& _xs,
        bool _intr
            )
        : beta(_beta), flmin(_flmin), X(_X), y(_y) 
        , q(_q), ju(_ju), vp(_vp), cl(_cl)
        , ne(_ne), nx(_nx), nlam(_nlam), maxit(_maxit), ulam(_ulam)
        , xb(_xb), xs(_xs)
        , intr(_intr)
        , g(_g)
        , a0(nlam)
        , ca(nx, nlam)
        , ia(nx)
        , kin(nlam)
        , dev(nlam)
        , alm(nlam)
    {
        a0.setZero();
        ca.setZero();
        ia.setZero();
        kin.setZero();
        dev.setZero();
        alm.setZero();
    }

    void fit() 
    {
        transl::spfishnet1<float, false>(
                beta, X, y, g, q, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, intr, maxit, xb, xs, lmu,
                a0, ca, ia, kin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() 
    {
        int ni = X.cols();
        int no = X.rows();
        int iisd = 0; // dummy: not used
        int iintr = intr;
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);

        spfishnet1_(const_cast<double*>(&beta), &no, &ni, 
                const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(),
                const_cast<double*>(y.data()),
                g.data(),
                const_cast<double*>(q.data()),
                const_cast<int*>(ju.data()), 
                const_cast<double*>(vp.data()), 
                const_cast<double*>(cl.data()),
                const_cast<int*>(&ne), 
                const_cast<int*>(&nx), 
                const_cast<int*>(&nlam),
                const_cast<double*>(&flmin), 
                const_cast<double*>(ulam.data()), 
                const_cast<double*>(&thr), 
                &iisd, &iintr,
                const_cast<int*>(&maxit), 
                const_cast<double*>(xb.data()),
                const_cast<double*>(xs.data()),
                &lmu,
                a0.data(), ca.data(), ia.data(), kin.data(), &dev0, dev.data(), alm.data(),
                &nlp, &jerr);
    }
};

struct sp_fishnet1_fixture
    : base_fixture
    , testing::WithParamInterface<
      std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, beta, flmin, isd, intr) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p);
        auto beta = dgen.make_beta(p);
        Eigen::VectorXd f = ((X * beta) / 2.).array().exp().matrix();
        std::mt19937 gen(seed);
        y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
        w = dgen.make_w(n);
        w.array() /= w.sum();
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        g.setZero(n);

        xm.resize(p);
        xs.resize(p);
        SpLStandardize2::eval(X, w, ju, isd, intr, xm, xs);
    }
protected:
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd cl;
    Eigen::VectorXd y, g, w, ulam, vp, xm, xs;
    Eigen::VectorXi ju;
    int nx, ne, maxit, nlam;
    bool isd, intr;
    double beta, flmin;

    void check_pack(const SpFishnet1Pack& actual,
                    const SpFishnet1Pack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        if (actual.jerr > util::max_active_reached_error().err_code(0)) {
            expect_near_vec(actual.g, expected.g, 3e-6);
        }
        expect_near_vec(actual.a0, expected.a0, 1e-8);
        expect_near_mat(actual.ca, expected.ca, 1e-7);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.kin, expected.kin);
        EXPECT_FLOAT_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.dev, expected.dev, 1e-7);
        expect_float_eq_vec(actual.alm, expected.alm);
    }
};

TEST_P(sp_fishnet1_fixture, sp_fishnet1_test)
{
    SpFishnet1Pack actual(
            beta, flmin, X, y, g, w, ju, vp, cl, ne, nx,
            nlam, maxit, ulam, xm, xs, intr);
    SpFishnet1Pack expected(actual);
    run(actual, expected, 10, 11);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpFishnet1Suite, sp_fishnet1_fixture,
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
