#include <legacy/legacy.h>
#include <testutil/translation/lognet.hpp>
#include <testutil/translation/splognet.hpp>
#include <translation/lognet_base_fixture.hpp>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct LognetDSCheckPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam, kopt;
    const double alpha, flmin;
    const bool isd, intr;

    Eigen::SparseMatrix<double> X;
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

    LognetDSCheckPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            int _kopt,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::MatrixXd& _g,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), kopt(_kopt), alpha(_alpha)
        , flmin(_flmin), isd(_isd), intr(_intr) 
        , X(_X)
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

    void fit() 
    {
        Eigen::MatrixXd X_d = X;
        transl::lognet<double, true>(
                alpha, X_d, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() 
    {
        // apply my fix and check if the results match with dense
        transl::splognet<double, true>(
                alpha, X, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }
};

struct lognet_dscheck_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, int> >
{
    void SetUp() override
    {
        size_t seed, n, p, nc;
        std::tie(seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt) = GetParam();
        DataGen dgen(seed); 
        X = dgen.make_X_sparse(n, p, 0.15);
        auto beta = dgen.make_beta(p);
        y.resize(n, std::max(static_cast<size_t>(2), nc));
        y.col(0) = dgen.make_y(X, beta);
        y.col(0).array() = (y.col(0).array() > y.col(0).mean()).template cast<double>();
        y.col(1).array() = 1. - y.col(0).array();
        jd = dgen.make_jd(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = p; // Fortran has an invalid read: must be p for testing.
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        g.setOnes(n, nc);
    }

protected:
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd cl, y, g;
    Eigen::VectorXd ulam, vp;
    Eigen::VectorXi jd;
    int nx, ne, maxit, nlam, kopt;
    double alpha, flmin;
    bool isd, intr; 

    void check_pack(const LognetDSCheckPack& actual,
                    const LognetDSCheckPack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);

        expect_eq_vec(actual.nin, expected.nin);
        expect_eq_vec(actual.ia, expected.ia);

        expect_near_vec(actual.dev, expected.dev, 1e-14);

        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            EXPECT_NEAR(actual.alm[i], expected.alm[i], actual.alm[i]*1e-15);
        }

        expect_near_mat(actual.g, expected.g, 4e-14);
        expect_near_mat(actual.y, expected.y, 1e-14);

        expect_near_mat(actual.a0, expected.a0, 1e-14);
        expect_near_vec(actual.ca, expected.ca, 1e-14);
    }
};

TEST_P(lognet_dscheck_fixture, lognet_dscheck_test)
{
    LognetDSCheckPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr, kopt,
            X, y, g, ulam, vp, cl, jd);
    LognetDSCheckPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    LognetDSCheckSuite, lognet_dscheck_fixture,
    testing::Combine(
        // seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt
        testing::Values(241, 3492, 3093772, 3827, 1993982874),
        testing::Values(10, 30, 50, 100),
        testing::Values(5, 20, 40),
        testing::Values(1, 2),
        testing::Values(50, 100),
        testing::Values(1, 10),
        testing::Values(1.0),
        testing::Values(1.0),
        testing::Bool(),
        testing::Bool(),
        testing::Values(0,1,2)
        )
);

} // namespace glmnetpp
