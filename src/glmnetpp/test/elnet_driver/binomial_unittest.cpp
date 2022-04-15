#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/thread.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/lognet.hpp>
#include <testutil/translation/splognet.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_driver/binomial.hpp>
#include <glmnetpp_bits/elnet_path.hpp>
#include <glmnetpp_bits/elnet_point.hpp>

namespace glmnetpp {

struct BinomialDriverPackBase
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam, kopt;
    const double alpha, flmin;
    const bool isd, intr;

    Eigen::MatrixXd y;
    Eigen::MatrixXd g;
    Eigen::MatrixXd cl;
    Eigen::VectorXd ca;
    Eigen::MatrixXd a0;
    Eigen::VectorXi ia;
    Eigen::VectorXi nin;
    Eigen::VectorXd dev; 
    Eigen::VectorXd alm;
    int nlp = 0, jerr = 0, lmu = 0;
    double dev0 = 0;
    
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXi& jd;

    BinomialDriverPackBase(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            int _kopt,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            const Eigen::MatrixXd& _y,
            const Eigen::MatrixXd& _g,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), kopt(_kopt), alpha(_alpha)
        , flmin(_flmin), isd(_isd), intr(_intr)
        , y(_y)
        , g(_g)
        , cl(_cl)
        , ca(_nx * _g.cols() * _nlam)
        , a0(_g.cols(), nlam)
        , ia(_nx)
        , nin(_nlam)
        , dev(_nlam)
        , alm(_nlam)
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

struct BinomialDriverPack : BinomialDriverPackBase
{
    using base_t = BinomialDriverPackBase;
    Eigen::MatrixXd X;
    BinomialDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            int _kopt, 
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::MatrixXd& _g,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _kopt, _alpha, _flmin, _isd, _intr, _y, _g, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::binomial>;
        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                alpha, X, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::lognet<double, true>(
                alpha, X, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }
};

struct binomial_driver_fixture_base
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, int> >
{
protected:
    static constexpr double tol = 1e-15;

    void check_pack(const BinomialDriverPackBase& actual,
                    const BinomialDriverPackBase& expected)
    {
        expect_near_mat(actual.g, expected.g, 8*tol);
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_mat(actual.a0, expected.a0, 2*tol);
        expect_near_vec(actual.ca, expected.ca, 2*tol);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.dev, expected.dev, 5*tol);

        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            EXPECT_NEAR(actual.alm[i], expected.alm[i], actual.alm[i]*1e-15);
        }

        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
        expect_near_mat(actual.y, expected.y, tol);
    }
};

struct binomial_driver_fixture : binomial_driver_fixture_base
{
protected:
    using base_t = binomial_driver_fixture_base;
    using base_t::tol;

    void check_pack(const BinomialDriverPack& actual,
                    const BinomialDriverPack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_mat(actual.X, expected.X);
    }
};

TEST_P(binomial_driver_fixture, binomial_driver_test)
{
    int seed, n, p, nc, maxit, nlam, kopt;
    double alpha, flmin;
    bool isd, intr;
    std::tie(seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X(n, p);
    auto beta = dgen.make_beta(p);
    Eigen::MatrixXd y(n, std::max(2, nc));
    y.col(0) = dgen.make_y(X, beta);
    y.col(0).array() = (y.col(0).array() > y.col(0).mean()).template cast<double>();
    y.col(1).array() = 1. - y.col(0).array();
    auto cl = dgen.make_cl(p);
    auto vp = dgen.make_vp(p);
    auto nx = dgen.make_nx(p);
    int ne = dgen.make_ne(p);
    auto ulam = dgen.make_ulam(nlam);
    auto jd = dgen.make_jd(p);
    Eigen::MatrixXd g(n, nc); g.setOnes();

    BinomialDriverPack actual(
            maxit, nx, ne, nlam, kopt, alpha, flmin, isd, intr,
            X, y, g, ulam, vp, cl, jd);
    BinomialDriverPack expected(actual);
    run(actual, expected, 14, 15); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    BinomialDriverSuite, binomial_driver_fixture,
    testing::Combine(
        // seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt
        testing::Values(241),
        testing::Values(10, 30, 50),
        testing::Values(5, 40, 60),
        testing::Values(1, 2),
        testing::Values(1, 50, 100),
        testing::Values(1, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool(),
        testing::Values(0, 1, 2)
        )
);

// =======================================================

struct SpBinomialDriverPack : BinomialDriverPackBase
{
    using base_t = BinomialDriverPackBase;
    Eigen::SparseMatrix<double> X;
    SpBinomialDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            int _kopt,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::MatrixXd& _g,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _kopt, _alpha, _flmin, _isd, _intr, _y, _g, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::binomial>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                alpha, X, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::splognet<double, true>(
                alpha, X, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }
};

struct sp_binomial_driver_fixture : binomial_driver_fixture_base
{
protected:
    using base_t = binomial_driver_fixture_base;
    using base_t::tol;

    void check_pack(const SpBinomialDriverPack& actual,
                    const SpBinomialDriverPack& expected)
    {
        Eigen::MatrixXd actual_X_dense = actual.X;
        Eigen::MatrixXd expected_X_dense = expected.X;
        expect_double_eq_mat(actual_X_dense, expected_X_dense);

        expect_near_mat(actual.g, expected.g, 5*tol);
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_mat(actual.a0, expected.a0, 2*tol);
        expect_near_vec(actual.ca, expected.ca, 2*tol);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.dev, expected.dev, 5*tol);

        // This check loosens expect_near_vec.
        // Do a relative check: if alm[i] is large, we want the relative error difference to be around 1e-15.
        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            EXPECT_NEAR(actual.alm[i], expected.alm[i], actual.alm[i]*1e-15);
        }

        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
    
        expect_near_mat(actual.y, expected.y, tol);
    }
};

TEST_P(sp_binomial_driver_fixture, sp_binomial_driver_test)
{
    int seed, n, p, nc, maxit, nlam, kopt;
    double alpha, flmin;
    bool isd, intr;
    std::tie(seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X_sparse(n, p);
    auto beta = dgen.make_beta(p);
    Eigen::MatrixXd y(n, std::max(2, nc));
    y.col(0) = dgen.make_y(X, beta);
    y.col(0).array() = (y.col(0).array() > y.col(0).mean()).template cast<double>();
    y.col(1).array() = 1. - y.col(0).array();
    auto jd = dgen.make_jd(p);
    auto vp = dgen.make_vp(p);
    auto cl = dgen.make_cl(p);
    auto nx = dgen.make_nx(p);
    int ne = dgen.make_ne(p);
    auto ulam = dgen.make_ulam(nlam);
    Eigen::MatrixXd g(n, nc); g.setOnes();

    SpBinomialDriverPack actual(
            maxit, nx, ne, nlam, kopt, alpha, flmin, isd, intr,
            X, y, g, ulam, vp, cl, jd);
    SpBinomialDriverPack expected(actual);
    run(actual, expected, 14, 15); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpBinomialDriverSuite, sp_binomial_driver_fixture,
    testing::Combine(
        // seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt
        testing::Values(241),
        testing::Values(10, 30, 50),
        testing::Values(5, 40, 60),
        testing::Values(1, 2),
        testing::Values(1, 50, 100),
        testing::Values(1, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool(),
        testing::Values(0, 1, 2)
        )
);

} // namespace glmnetpp
