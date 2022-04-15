#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/translation/fishnet1.hpp>
#include <testutil/thread.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <glmnetpp_bits/elnet_path/poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/poisson_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>

namespace glmnetpp {

struct PoissonNaivePack
{
    const double beta, flmin, thr = 1e-14;
    const Eigen::MatrixXd& X;
    const Eigen::VectorXd& y;
    const Eigen::VectorXd& q;
    const Eigen::VectorXi& ju;
    const Eigen::VectorXd& vp;
    const Eigen::MatrixXd& cl;
    const int ne, nx, nlam, maxit;
    const Eigen::VectorXd& ulam;
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

    PoissonNaivePack(
        double _beta, 
        double _flmin, 
        const Eigen::MatrixXd& _X,
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
        bool _intr
            )
        : beta(_beta), flmin(_flmin), X(_X), y(_y) 
        , q(_q), ju(_ju), vp(_vp), cl(_cl)
        , ne(_ne), nx(_nx), nlam(_nlam), maxit(_maxit), ulam(_ulam)
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
        using internal_t = ElnetPointInternal<
                    util::glm_type::poisson,
                    util::mode_type<util::glm_type::poisson>::naive,
                    double,
                    int, 
                    int>;
        using elnet_point_t = ElnetPoint<
                util::glm_type::poisson,
                util::mode_type<util::glm_type::poisson>::naive,
                internal_t>;
        using elnet_path_t = ElnetPath<
            util::glm_type::poisson,
            util::mode_type<util::glm_type::poisson>::naive,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                beta, ju, vp, cl, ne, nx, X, y, g, q, 
                nlam, flmin, ulam, thr, intr, maxit, lmu,
                a0, ca, ia, kin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() 
    {
        transl::fishnet1<double>(
                beta, X, y, g, q, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, intr, maxit, lmu,
                a0, ca, ia, kin, dev0, dev, alm, nlp, jerr);
    }
};

struct poisson_naive_fixture
    : base_fixture
    , testing::WithParamInterface<
      std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        bool do_offset;
        std::tie(seed, n, p, maxit, nlam, beta, flmin, do_offset, isd, intr) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X(n, p);
        auto beta = dgen.make_beta(p, 0.8);
        w = dgen.make_w(n);
        w.array() /= w.sum();
        Eigen::VectorXd f = dgen.make_y(X, beta);
        double fm = f.dot(w);
        double fs = std::sqrt((f.array() - fm).square().matrix().dot(w));
        f.array() -= fm;
        f /= fs;
        std::mt19937 gen(seed);
        y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        if (do_offset) {
            g = g.NullaryExpr(n, [&](auto) { return 0.5 + std::normal_distribution<double>()(gen); });
        } else {
            g.setZero(n);
        }

        Eigen::VectorXd xm(p), xs(p);
        LStandardize1::eval(X, w, ju, isd, intr, xm, xs);
    }
protected:
    Eigen::MatrixXd X, cl;
    Eigen::VectorXd y, g, w, ulam, vp;
    Eigen::VectorXi ju;
    int nx, ne, maxit, nlam;
    bool isd, intr;
    double beta, flmin;

    void check_pack(const PoissonNaivePack& actual,
                    const PoissonNaivePack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);

        // only check offset if this error did not occur
        // Fortran changes the order of the active set update which does an extra update of offset.
        if (actual.jerr > util::max_active_reached_error().err_code(0)) {
            expect_near_vec(actual.g, expected.g, actual.g.array().abs().maxCoeff() * 3e-15);
        }

        expect_near_vec(actual.a0, expected.a0, actual.a0.array().abs().maxCoeff() * 2e-15);
        expect_near_mat(actual.ca, expected.ca, 2e-15);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.kin, expected.kin);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.dev, expected.dev, 2e-15);
        expect_double_eq_vec(actual.alm, expected.alm);
    }
};

TEST_P(poisson_naive_fixture, poisson_naive_test)
{
    PoissonNaivePack actual(
            beta, flmin, X, y, g, w, ju, vp, cl, ne, nx,
            nlam, maxit, ulam, intr);
    PoissonNaivePack expected(actual);
    run(actual, expected, 10, 11);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    PoissonNaiveSuite, poisson_naive_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0),            // flmin
        testing::Bool(),                      // make offset or not
        testing::Bool(),                      // isd
        testing::Bool()                       // intr
        )
);

} // namespace glmnetpp
