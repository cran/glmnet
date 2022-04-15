#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/thread.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/fishnet.hpp>
#include <testutil/translation/spfishnet.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_driver/poisson.hpp>
#include <glmnetpp_bits/elnet_path.hpp>
#include <glmnetpp_bits/elnet_point.hpp>

namespace glmnetpp {

struct PoissonDriverPackBase
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, intr;

    Eigen::VectorXd g;
    Eigen::MatrixXd cl;
    Eigen::MatrixXd ca;
    Eigen::VectorXd a0;
    Eigen::VectorXi ia;
    Eigen::VectorXi nin;
    Eigen::VectorXd dev; 
    Eigen::VectorXd alm;
    int nlp = 0, jerr = 0, lmu = 0;
    double dev0 = 0;
    
    const Eigen::VectorXd& y;
    const Eigen::VectorXd& w;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXi& jd;

    PoissonDriverPackBase(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha)
        , flmin(_flmin), isd(_isd), intr(_intr) 
        , g(_g)
        , cl(_cl)
        , ca(_nx, _nlam)
        , a0(nlam)
        , ia(_nx)
        , nin(_nlam)
        , dev(_nlam)
        , alm(_nlam)
        , y(_y)
        , w(_w)
        , ulam(_ulam)
        , vp(_vp)
        , jd(_jd)
    {
        ca.setZero();
        a0.setZero();
        ia.setZero();
        nin.setZero();
        dev.setZero();
        alm.setZero();
    }

};

struct PoissonDriverPack : PoissonDriverPackBase
{
    using base_t = PoissonDriverPackBase;
    Eigen::MatrixXd X;
    PoissonDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            const Eigen::MatrixXd& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _isd, _intr, _y, _g, _w, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::poisson>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                alpha, X, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::fishnet<double>(
                alpha, X, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }
};

struct poisson_driver_fixture_base
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool> >
{
protected:
    static constexpr double tol = 1e-15;

    void check_pack(const PoissonDriverPackBase& actual,
                    const PoissonDriverPackBase& expected)
    {
        // only check offset if this error did not occur
        // Fortran changes the order of the active set update which does an extra update of offset.
        if (actual.jerr > util::max_active_reached_error().err_code(0)) {
            expect_near_vec(actual.g, expected.g, 8*tol);
        }
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_vec(actual.a0, expected.a0, 3*tol);
        expect_near_mat(actual.ca, expected.ca, 3*tol);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.dev, expected.dev, 5*tol);
        expect_double_eq_vec(actual.alm, expected.alm);
        EXPECT_NEAR(actual.nlp, expected.nlp, 1); // nlp might be off by 1 if y's are large
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.y, expected.y, tol);
    }
};

struct poisson_driver_fixture : poisson_driver_fixture_base
{
protected:
    using base_t = poisson_driver_fixture_base;
    using base_t::tol;

    void check_pack(const PoissonDriverPack& actual,
                    const PoissonDriverPack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_mat(actual.X, expected.X);
    }
};

TEST_P(poisson_driver_fixture, poisson_driver_test)
{
    int seed, n, p, maxit, nlam;
    double alpha, flmin;
    bool isd, intr;
    std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr) = GetParam();
    DataGen dgen(seed); 
    auto X = dgen.make_X(n, p);
    auto beta = dgen.make_beta(p);
    Eigen::VectorXd f = ((X * beta) / 2.).array().exp().matrix();
    std::mt19937 gen(seed);
    Eigen::VectorXd y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
    auto w = dgen.make_w(n);
    auto jd = dgen.make_jd(p);
    auto vp = dgen.make_vp(p);
    auto cl = dgen.make_cl(p);
    auto nx = dgen.make_nx(p);
    auto ne = dgen.make_ne(p);
    auto ulam = dgen.make_ulam(nlam);
    Eigen::VectorXd g; g.setZero(n);

    PoissonDriverPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr,
            X, y, g, w, ulam, vp, cl, jd);
    PoissonDriverPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    PoissonDriverSuite, poisson_driver_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr
        testing::Values(241, 132, 5241, 231),
        testing::Values(10, 30, 50),
        testing::Values(5, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool()
        )
);

// =======================================================

struct SpPoissonDriverPack : PoissonDriverPackBase
{
    using base_t = PoissonDriverPackBase;
    const Eigen::SparseMatrix<double>& X;
    SpPoissonDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _isd, _intr, _y, _g, _w, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::poisson>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                alpha, X, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::spfishnet<double, true>(
                alpha, X, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }
};

struct sp_poisson_driver_fixture : poisson_driver_fixture_base
{
protected:
    using base_t = poisson_driver_fixture_base;
    using base_t::tol;

    void check_pack(const SpPoissonDriverPack& actual,
                    const SpPoissonDriverPack& expected)
    {
        expect_near_vec(actual.g, expected.g, 11*tol);
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_vec(actual.a0, expected.a0, 6*tol);
        expect_near_mat(actual.ca, expected.ca, 4*tol);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.dev, expected.dev, 5*tol);

        // This check loosens expect_near_vec.
        // When alm is large (>= 1), check for exact equality.
        // Otherwise, if too small, put a tolerance of 1e-15.
        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            if (actual.alm[i] < 1) {
                EXPECT_NEAR(actual.alm[i], expected.alm[i], 1e-15);
            } else {
                EXPECT_DOUBLE_EQ(actual.alm[i], expected.alm[i]);
            }
        }

        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
    
        expect_near_vec(actual.y, expected.y, tol);
    }
};

TEST_P(sp_poisson_driver_fixture, sp_poisson_driver_test)
{
    int seed, n, p, maxit, nlam;
    double alpha, flmin;
    bool isd, intr;
    std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr) = GetParam();
    DataGen dgen(seed); 
    auto X = dgen.make_X_sparse(n, p);
    auto beta = dgen.make_beta(p, 0.2);
    Eigen::VectorXd f = ((X * beta) / 2.).array().exp().matrix();
    std::mt19937 gen(seed);
    Eigen::VectorXd y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
    auto w = dgen.make_w(n);
    auto jd = dgen.make_jd(p);
    auto vp = dgen.make_vp(p);
    auto cl = dgen.make_cl(p);
    auto nx = dgen.make_nx(p);
    auto ne = dgen.make_ne(p);
    auto ulam = dgen.make_ulam(nlam);
    Eigen::VectorXd g(n); g.setZero();

    SpPoissonDriverPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr,
            X, y, g, w, ulam, vp, cl, jd);
    SpPoissonDriverPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpPoissonDriverSuite, sp_poisson_driver_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr
        testing::Values(241, 132, 5241, 231),
        testing::Values(10, 30, 50),
        testing::Values(5, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool()
        )
);

} // namespace glmnetpp
