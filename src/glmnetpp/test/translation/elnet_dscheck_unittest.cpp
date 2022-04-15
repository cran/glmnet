#include <legacy/legacy.h>
#include <testutil/translation/elnet.hpp>
#include <testutil/translation/spelnet.hpp>
#include <translation/elnet_base_fixture.hpp>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct ElnetDSCheckPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, intr, ka;

    Eigen::SparseMatrix<double> X;
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

    ElnetDSCheckPack(
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
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha)
        , flmin(_flmin), isd(_isd), intr(_intr), ka(_ka) 
        , X(_X)
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

    void fit() 
    {
        Eigen::MatrixXd X_d = X;
        transl::elnet<double>(
                ka, alpha, X_d, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr);
    }

    void fit_old() 
    {
        // apply my fix and check if the results match with dense
        transl::spelnet<double>(
                ka, alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr);
    }
};

struct elnet_dscheck_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p);
        auto beta = dgen.make_beta(p);
        y = dgen.make_y(X, beta);
        w = dgen.make_w(n);
        jd = dgen.make_jd(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
    }

protected:
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd cl;
    Eigen::VectorXd y, w, ulam, vp;
    Eigen::VectorXi jd;
    int nx, ne, maxit, nlam;
    double alpha, flmin;
    bool isd, intr, ka;

    void check_pack(const ElnetDSCheckPack& actual,
                    const ElnetDSCheckPack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        expect_near_vec(actual.rsq, expected.rsq, 1e-14);

        expect_eq_vec(actual.nin, expected.nin);
        expect_eq_vec(actual.ia, expected.ia);

        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            EXPECT_NEAR(actual.alm[i], expected.alm[i], actual.alm[i]*1e-15);
        }

        expect_near_vec(actual.a0, expected.a0, 1e-14);
        expect_near_mat(actual.ca, expected.ca, 1e-14);

        // Note: we do not compare y because dense and sparse treat it differently.
        // In dense, it is the full residual. 
        // In sparse, it's only a partial residual coming from the raw (uncentered) x_i
        // (i.e. it is missing the component that comes from centering of x_i).
    }
};

TEST_P(elnet_dscheck_fixture, elnet_dscheck_test)
{
    ElnetDSCheckPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr, ka,
            X, y, w, ulam, vp, cl, jd);
    ElnetDSCheckPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    ElnetDSCheckSuite, elnet_dscheck_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka
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
