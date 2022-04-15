#include <legacy/legacy.h>
#include <testutil/translation/lognet.hpp>
#include <translation/lognet_base_fixture.hpp>

namespace glmnetpp {

struct LognetPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam, kopt;
    const double alpha, flmin;
    const bool isd, intr;

    Eigen::MatrixXd X;
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

    LognetPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            bool _isd,
            bool _intr,
            int _kopt,
            const Eigen::MatrixXd& _X,
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
        transl::lognet<float, false>(
                alpha, X, y, g, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old()
    {
        int no = X.rows();
        int ni = X.cols();
        int nc = g.cols();
        int iisd = isd;
        int iintr = intr;
        ::lognet_(
                const_cast<double*>(&alpha), &no, &ni, &nc,
                X.data(), y.data(), g.data(),
                const_cast<int*>(jd.data()),
                const_cast<double*>(vp.data()),
                cl.data(),
                const_cast<int*>(&ne), 
                const_cast<int*>(&nx), 
                const_cast<int*>(&nlam), 
                const_cast<double*>(&flmin), 
                const_cast<double*>(ulam.data()), 
                const_cast<double*>(&thr), &iisd, &iintr, 
                const_cast<int*>(&maxit), 
                const_cast<int*>(&kopt), &lmu, a0.data(), ca.data(), ia.data(), 
                nin.data(), &dev0, dev.data(), alm.data(), &nlp, &jerr);
    }
};

struct lognet_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, int> >
{
    void SetUp() override
    {
        size_t seed, n, p, nc;
        std::tie(seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt) = GetParam();
        DataGen dgen(seed); 
        X = dgen.make_X(n, p);
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
    Eigen::MatrixXd X, cl, y, g;
    Eigen::VectorXd ulam, vp;
    Eigen::VectorXi jd;
    int nx, ne, maxit, nlam, kopt;
    double alpha, flmin;
    bool isd, intr; 

    void check_pack(const LognetPack& actual,
                    const LognetPack& expected)
    {
        expect_float_eq_mat(actual.X, expected.X);
        expect_near_mat(actual.g, expected.g, 1e-7);

        // Fortran version uses uninitialized values in columns that are omitted (inside ju).
        // My version sets them to 0.
        // As a heuristic, if my cl contains (0,0) in a column, 
        // we will assume that's an omitted feature, so we skip the check.
        // Otherwise, valgrind will complain.
        for (int j = 0; j < actual.cl.cols(); ++j) {
            if ((actual.cl.col(j).array() == 0).all()) continue;
            for (int i = 0; i < actual.cl.rows(); ++i) {
                EXPECT_FLOAT_EQ(actual.cl(i,j), expected.cl(i,j));
            }
        }

        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.nin, expected.nin);

        if (actual.a0.cols() > 0) {
            expect_near_mat(actual.a0.block(actual.a0.rows(), actual.a0.cols()-1, 0, 0), 
                            expected.a0.block(expected.a0.rows(), expected.a0.cols()-1, 0, 0), 1e-7);
        }
        int nx = actual.nx;
        int nc = actual.g.cols();
        int nlam = actual.nlam;
        for (int i = 0; i < nlam - 1; ++i) {
            Eigen::Map<const Eigen::MatrixXd> actual_a_slice(
                    actual.ca.data() + nx * nc * i, nx, nc);
            Eigen::Map<const Eigen::MatrixXd> expected_a_slice(
                    expected.ca.data() + nx * nc * i, nx, nc);
            expect_near_mat(actual_a_slice, expected_a_slice, 1e-7);
        }

        expect_float_eq_mat(actual.y, expected.y);
        EXPECT_FLOAT_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.dev, expected.dev, 1e-7);
        expect_float_eq_vec(actual.alm, expected.alm);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);

    }
};

TEST_P(lognet_fixture, lognet_test)
{
    LognetPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr, kopt,
            X, y, g, ulam, vp, cl, jd);
    LognetPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    LognetSuite, lognet_fixture,
    testing::Combine(
        // seed, n, p, nc, maxit, nlam, alpha, flmin, isd, intr, kopt
        testing::Values(241),
        testing::Values(10, 30, 50),
        testing::Values(5, 20, 40),
        testing::Values(1, 2),
        testing::Values(1, 50, 100),
        testing::Values(1, 4),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool(),
        testing::Values(0, 1, 2)
        )
);

} // namespace glmnetpp
