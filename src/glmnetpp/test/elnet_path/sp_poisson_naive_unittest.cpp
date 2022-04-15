#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/translation/spfishnet1.hpp>
#include <testutil/thread.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <glmnetpp_bits/elnet_path/sp_poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/sp_poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp> // naive
#include <Eigen/SparseCore>

namespace glmnetpp {

struct SpPoissonNaivePack
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

    SpPoissonNaivePack(
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
        using internal_t = SpElnetPointInternal<
                    util::glm_type::poisson,
                    util::mode_type<util::glm_type::poisson>::naive,
                    double,
                    int, 
                    int>;
        using elnet_point_t = SpElnetPoint<
                util::glm_type::poisson,
                util::mode_type<util::glm_type::poisson>::naive,
                internal_t>;
        using elnet_path_t = SpElnetPath<
            util::glm_type::poisson,
            util::mode_type<util::glm_type::poisson>::naive,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                beta, ju, vp, cl, ne, nx, X, y, g, q, 
                nlam, flmin, ulam, xb, xs, thr, intr, maxit, lmu,
                a0, ca, ia, kin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() 
    {
        transl::spfishnet1<double, true>(
                beta, X, y, g, q, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, intr, maxit, xb, xs, lmu,
                a0, ca, ia, kin, dev0, dev, alm, nlp, jerr);
    }
};

struct sp_poisson_naive_fixture
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

    void check_pack(const SpPoissonNaivePack& actual,
                    const SpPoissonNaivePack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        expect_near_vec(actual.g, expected.g, 4e-15);
        expect_near_vec(actual.a0, expected.a0, 6e-15);
        expect_near_mat(actual.ca, expected.ca, 2e-15);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.kin, expected.kin);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.dev, expected.dev, 2e-15);
        expect_double_eq_vec(actual.alm, expected.alm);
    }
};

TEST_P(sp_poisson_naive_fixture, sp_poisson_naive_test)
{
    SpPoissonNaivePack actual(
            beta, flmin, X, y, g, w, ju, vp, cl, ne, nx,
            nlam, maxit, ulam, xm, xs, intr);
    SpPoissonNaivePack expected(actual);
    run(actual, expected, 10, 11);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpPoissonNaiveSuite, sp_poisson_naive_fixture,
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
