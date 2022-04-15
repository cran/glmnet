#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/thread.hpp>
#include <legacy/legacy.h>
#include <testutil/translation/fishnet.hpp>

namespace glmnetpp {

struct FishnetPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, intr;

    Eigen::MatrixXd X;
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

    FishnetPack(
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
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha)
        , flmin(_flmin), isd(_isd), intr(_intr) 
        , X(_X)
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

    void fit()
    {
        transl::fishnet<float>(
                alpha, X, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old()
    {
        int no = X.rows();
        int ni = X.cols();
        int iisd = isd;
        int iintr = intr;
        ::fishnet_(
                const_cast<double*>(&alpha), &no, &ni,
                X.data(), 
                const_cast<double*>(y.data()), g.data(), 
                const_cast<double*>(w.data()),
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
                &lmu, a0.data(), ca.data(), ia.data(), 
                nin.data(), &dev0, dev.data(), alm.data(), &nlp, &jerr);
    }
};

struct fishnet_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr) = GetParam();
        DataGen dgen(seed); 
        X = dgen.make_X(n, p);
        auto beta = dgen.make_beta(p);
        Eigen::VectorXd f = X * beta;
        f.array() /= std::max(beta.norm(), 1e-4);
        f = f.array().exp().matrix();
        std::mt19937 gen(seed);
        y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
        w = dgen.make_w(n);
        jd = dgen.make_jd(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        g.setZero(n);
    }

protected:
    Eigen::MatrixXd X, cl;
    Eigen::VectorXd y, g, w, ulam, vp;
    Eigen::VectorXi jd;
    int nx, ne, maxit, nlam;
    double alpha, flmin;
    bool isd, intr; 

    void check_pack(const FishnetPack& actual,
                    const FishnetPack& expected)
    {
        expect_float_eq_mat(actual.X, expected.X);
        expect_near_vec(actual.g, expected.g, 1e-7);

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

        expect_near_vec(actual.a0, expected.a0, 1e-7);
        expect_near_mat(actual.ca, expected.ca, 1e-7);

        expect_float_eq_mat(actual.y, expected.y);
        EXPECT_FLOAT_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.dev, expected.dev, 1e-7);
        expect_float_eq_vec(actual.alm, expected.alm);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
    }
};

TEST_P(fishnet_fixture, fishnet_test)
{
    FishnetPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr,
            X, y, g, w, ulam, vp, cl, jd);
    FishnetPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    FishnetSuite, fishnet_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr
        testing::Values(241),
        testing::Values(10, 30, 50),
        testing::Values(5, 20, 40),
        testing::Values(1, 50, 100),
        testing::Values(1, 4),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool()
        )
);

} // namespace glmnetpp
