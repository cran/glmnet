#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/thread.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/elnet.hpp>
#include <testutil/translation/spelnet.hpp>
#include <testutil/translation/multelnet.hpp>
#include <testutil/translation/multspelnet.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_driver/gaussian.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_multi.hpp>

namespace glmnetpp {

struct GaussianDriverPackBase
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, intr, ka;

    Eigen::VectorXd y;
    Eigen::VectorXd w;
    Eigen::MatrixXd cl;
    Eigen::MatrixXd ca;
    Eigen::VectorXd a0;
    Eigen::VectorXi ia;
    Eigen::VectorXi nin;
    Eigen::VectorXd rsq; 
    Eigen::VectorXd alm;
    int nlp = 0, jerr = 0, lmu = 0;
    
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXi& jd;

    GaussianDriverPackBase(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            bool _ka,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha)
        , flmin(_flmin), isd(_isd), intr(_intr), ka(_ka) 
        , y(_y)
        , w(_w)
        , cl(_cl)
        , ca(_nx, _nlam)
        , a0(nlam)
        , ia(_nx)
        , nin(_nlam)
        , rsq(_nlam)
        , alm(_nlam)
        , ulam(_ulam)
        , vp(_vp)
        , jd(_jd)
    {
        ca.setZero();
        a0.setZero();
        ia.setZero();
        nin.setZero();
        rsq.setZero();
        alm.setZero();
    }

};

struct GaussianDriverPack : GaussianDriverPackBase
{
    using base_t = GaussianDriverPackBase;
    Eigen::MatrixXd X;
    GaussianDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            bool _ka,
            const Eigen::MatrixXd& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _isd, _intr, _ka, _y, _w, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::gaussian>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                ka, alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::elnet<double>(
                ka, alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr);
    }
};

struct gaussian_driver_fixture_base
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
protected:
    static constexpr double tol = 1e-15;

    void check_pack(const GaussianDriverPackBase& actual,
                    const GaussianDriverPackBase& expected)
    {
        expect_double_eq_vec(actual.w, expected.w);
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_vec(actual.a0, expected.a0, tol);
        expect_near_mat(actual.ca, expected.ca, tol);
        expect_near_vec(actual.ia, expected.ia, 1);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.rsq, expected.rsq, tol);
        expect_double_eq_vec(actual.alm, expected.alm);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
    
        expect_near_vec(actual.y, expected.y, tol);
    }
};

struct gaussian_driver_fixture : gaussian_driver_fixture_base
{
protected:
    using base_t = gaussian_driver_fixture_base;
    using base_t::tol;

    void check_pack(const GaussianDriverPack& actual,
                    const GaussianDriverPack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_mat(actual.X, expected.X);
    }
};

TEST_P(gaussian_driver_fixture, gaussian_driver_test)
{
    int seed, n, p, maxit, nlam;
    double alpha, flmin;
    bool isd, intr, ka;
    std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X(n, p);
    auto beta = dgen.make_beta(p);
    auto y = dgen.make_y(X, beta);
    auto w = dgen.make_w(n);
    auto cl = dgen.make_cl(p);
    auto vp = dgen.make_vp(p);
    auto nx = dgen.make_nx(p);
    int ne = dgen.make_ne(p);
    auto ulam = dgen.make_ulam(nlam);
    auto jd = dgen.make_jd(p);

    GaussianDriverPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr, ka,
            X, y, w, ulam, vp, cl, jd);
    GaussianDriverPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    GaussianDriverSuite, gaussian_driver_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka
        testing::Values(241, 412, 23968),
        testing::Values(10, 30, 50),
        testing::Values(5, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

// =======================================================

struct SpGaussianDriverPack : GaussianDriverPackBase
{
    using base_t = GaussianDriverPackBase;
    Eigen::SparseMatrix<double> X;
    SpGaussianDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            bool _ka,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _isd, _intr, _ka, _y, _w, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::gaussian>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                ka, alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::spelnet<double>(
                ka, alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr);
    }
};

struct sp_gaussian_driver_fixture : gaussian_driver_fixture_base
{
protected:
    using base_t = gaussian_driver_fixture_base;
    using base_t::tol;

    void check_pack(const SpGaussianDriverPack& actual,
                    const SpGaussianDriverPack& expected)
    {
        Eigen::MatrixXd actual_X_dense = actual.X;
        Eigen::MatrixXd expected_X_dense = expected.X;
        expect_double_eq_mat(actual_X_dense, expected_X_dense);

        expect_double_eq_vec(actual.w, expected.w);
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_vec(actual.a0, expected.a0, tol);
        expect_near_mat(actual.ca, expected.ca, 1e-14);
        expect_near_vec(actual.ia, expected.ia, 1);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.rsq, expected.rsq, 1e-14);

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
    
        expect_near_vec(actual.y, expected.y, 1e-14);
    }
};

TEST_P(sp_gaussian_driver_fixture, sp_gaussian_driver_test)
{
    int seed, n, p, maxit, nlam;
    double alpha, flmin;
    bool isd, intr, ka;
    std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X_sparse(n, p);
    auto beta = dgen.make_beta(p);
    auto y = dgen.make_y(X, beta);
    auto w = dgen.make_w(n);
    auto jd = dgen.make_jd(p);
    auto vp = dgen.make_vp(p);
    auto cl = dgen.make_cl(p);
    auto nx = dgen.make_nx(p);
    int ne = dgen.make_ne(p);
    auto ulam = dgen.make_ulam(nlam);

    SpGaussianDriverPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr, ka,
            X, y, w, ulam, vp, cl, jd);
    SpGaussianDriverPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpGaussianDriverSuite, sp_gaussian_driver_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka
        testing::Values(241, 412, 23968),
        testing::Values(10, 30, 50),
        testing::Values(5, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

// =======================================================
// Multi-response Driver
// =======================================================

struct GaussianMultiDriverPackBase
    : GaussianDriverPackBase
{
    using base_t = GaussianDriverPackBase;

    const bool jsd;

    // will be modified
    Eigen::MatrixXd y;
    Eigen::MatrixXd a0;
    Eigen::VectorXd ao;

    GaussianMultiDriverPackBase(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            bool _isd,
            bool _jsd,
            bool _intr,
            double _alpha,
            double _flmin,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _isd, _intr, false /* not used */, 
                _y.col(0) /* not used */, _w, _ulam, _vp, _cl, _jd)
        , jsd(_jsd)
        , y(_y)
        , a0(_y.cols(), _nlam)
        , ao(_nx * _y.cols() * _nlam)
    {
        a0.setZero();
        ao.setZero();
    }
};

struct GaussianMultiDriverPack : GaussianMultiDriverPackBase
{
    using base_t = GaussianMultiDriverPackBase;
    Eigen::MatrixXd X;
    GaussianMultiDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            bool _isd,
            bool _jsd,
            bool _intr,
            double _alpha,
            double _flmin,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _isd, _jsd, _intr, _alpha, _flmin, _y, _w, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::gaussian>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                alpha, X, y, w, jd, vp, cl, ne, nx, 
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::multelnet<double>(
                alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, nin, rsq, alm, nlp, jerr);
    }
};

struct gaussian_multi_driver_fixture_base
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
protected:
    static constexpr double tol = 1e-15;

    void check_pack(const GaussianMultiDriverPackBase& actual,
                    const GaussianMultiDriverPackBase& expected)
    {
        expect_double_eq_vec(actual.w, expected.w);
        expect_near_mat(actual.a0, expected.a0, tol);
        expect_near_vec(actual.ao, expected.ao, tol);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.rsq, expected.rsq, tol);
        expect_double_eq_vec(actual.alm, expected.alm);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
        expect_near_mat(actual.y, expected.y, tol);
    }
};

struct gaussian_multi_driver_fixture : gaussian_multi_driver_fixture_base
{
    void SetUp() override
    {
        size_t seed, n, p, nr = 3;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, jsd, intr) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X(n, p);
        y.resize(n, nr);
        for (size_t i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        w = dgen.make_w(n);
        jd = dgen.make_jd(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
    }

protected:
    using base_t = gaussian_multi_driver_fixture_base;
    using base_t::tol;
    Eigen::MatrixXd X, y, cl;
    Eigen::VectorXd w, vp, ulam;
    Eigen::VectorXi jd;
    size_t maxit, nlam, nx, ne;
    double alpha, flmin;
    bool isd, jsd, intr;

    void check_pack(const GaussianMultiDriverPack& actual,
                    const GaussianMultiDriverPack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_mat(actual.X, expected.X);
    }
};

TEST_P(gaussian_multi_driver_fixture, gaussian_multi_driver_test)
{
    GaussianMultiDriverPack actual(
            maxit, nx, ne, nlam, isd, jsd, intr, alpha, flmin,
            X, y, w, ulam, vp, cl, jd);
    GaussianMultiDriverPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    GaussianMultiDriverSuite, gaussian_multi_driver_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5),       // flmin
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

// =======================================================

struct SpGaussianMultiDriverPack : GaussianMultiDriverPackBase
{
    using base_t = GaussianMultiDriverPackBase;
    const Eigen::SparseMatrix<double> X;
    SpGaussianMultiDriverPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            bool _isd,
            bool _jsd,
            bool _intr,
            double _alpha,
            double _flmin,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : base_t(_maxit, _nx, _ne, _nlam, _isd, _jsd, _intr, _alpha, _flmin, _y, _w, _ulam, _vp, _cl, _jd)
        , X(_X)
    {}

    void fit()
    {
        using elnet_driver_t = ElnetDriver<
            util::glm_type::gaussian>;

        elnet_driver_t elnet_driver;
        elnet_driver.fit(
                alpha, X, y, w, jd, vp, cl, ne, nx, 
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old()
    {
        transl::multspelnet<double>(
                alpha, X, y, w, jd, vp, cl, ne, nx, 
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, nin, rsq, alm, nlp, jerr);
    }
};

struct sp_gaussian_multi_driver_fixture : gaussian_multi_driver_fixture_base
{
    void SetUp() override
    {
        size_t seed, n, p, nr = 3;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, jsd, intr) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p);
        y.resize(n, nr);
        for (size_t i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        w = dgen.make_w(n);
        jd = dgen.make_jd(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
    }

protected:
    using base_t = gaussian_multi_driver_fixture_base;
    using base_t::tol;
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd y, cl;
    Eigen::VectorXd w, vp, ulam;
    Eigen::VectorXi jd;
    size_t maxit, nlam, nx, ne;
    double alpha, flmin;
    bool isd, jsd, intr;

    void check_pack(const SpGaussianMultiDriverPack& actual,
                    const SpGaussianMultiDriverPack& expected)
    {

        expect_double_eq_vec(actual.w, expected.w);
        expect_double_eq_mat(actual.cl, expected.cl);
        expect_near_mat(actual.a0, expected.a0, tol);
        expect_near_vec(actual.ao, expected.ao, 1e-15);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);
        expect_near_vec(actual.rsq, expected.rsq, 2e-15);

        // This check loosens expect_near_vec.
        // The absolute difference proportion to the actual value must be within a certain tolerance.
        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            EXPECT_NEAR(actual.alm[i], expected.alm[i], 
                    actual.alm[i] * 1e-15);
        }

        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
    
        expect_near_mat(actual.y, expected.y, 4e-15);
    }
};

TEST_P(sp_gaussian_multi_driver_fixture, sp_gaussian_multi_driver_test)
{
    SpGaussianMultiDriverPack actual(
            maxit, nx, ne, nlam, isd, jsd, intr, alpha, flmin,
            X, y, w, ulam, vp, cl, jd);
    SpGaussianMultiDriverPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpGaussianMultiDriverSuite, sp_gaussian_multi_driver_fixture,
    testing::Combine(
        testing::Values(241, 412, 3, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5),       // flmin
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);


} // namespace glmnetpp
