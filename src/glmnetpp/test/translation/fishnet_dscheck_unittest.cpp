#include <legacy/legacy.h>
#include <Eigen/SparseCore>
#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <testutil/translation/fishnet.hpp>
#include <testutil/translation/spfishnet.hpp>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct FishnetDSCheckPack
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
    
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& y;
    const Eigen::VectorXd& w;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXi& jd;

    FishnetDSCheckPack(
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
        , X(_X)
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
        Eigen::MatrixXd X_d = X;
        transl::fishnet<double>(
                alpha, X_d, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() 
    {
        // apply my fix and check if the results match with dense
        transl::spfishnet<double, true>(
                alpha, X, y, g, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, dev0, dev, alm, nlp, jerr);
    }
};

struct fishnet_dscheck_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        bool do_offset;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr, do_offset) = GetParam();
        DataGen dgen(seed); 
        X = dgen.make_X_sparse(n, p);
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
        nx = p;
        ne = p+1;
        ulam = dgen.make_ulam(nlam);
        if (do_offset) {
            g.setOnes(n);
        } else {
            g.setZero(n);
        }
    }

protected:
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd cl;
    Eigen::VectorXd y, g, w, ulam, vp;
    Eigen::VectorXi jd;
    int nx, ne, maxit, nlam;
    double alpha, flmin;
    bool isd, intr; 

    void check_pack(const FishnetDSCheckPack& actual,
                    const FishnetDSCheckPack& expected)
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

        expect_near_vec(actual.g, expected.g, 4e-14);
        expect_near_vec(actual.a0, expected.a0, 1e-14);
        expect_near_mat(actual.ca, expected.ca, 1e-14);
    }
};

TEST_P(fishnet_dscheck_fixture, fishnet_dscheck_test)
{
    FishnetDSCheckPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr,
            X, y, g, w, ulam, vp, cl, jd);
    FishnetDSCheckPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    FishnetDSCheckSuite, fishnet_dscheck_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr, do_offset
        testing::Values(241, 3492, 3093772, 3827, 1993982874),
        testing::Values(10, 30, 50, 100),
        testing::Values(5, 20, 40),
        testing::Values(50, 100),
        testing::Values(1, 10),
        testing::Values(1.0),
        testing::Values(1.0),
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

} // namespace glmnetpp
