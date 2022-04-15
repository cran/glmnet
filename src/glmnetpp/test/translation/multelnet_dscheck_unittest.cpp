#include <legacy/legacy.h>
#include <testutil/translation/multelnet.hpp>
#include <testutil/translation/multspelnet.hpp>
#include <translation/elnet_base_fixture.hpp>

// separately unit-tested
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct MultElnetDSCheckPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, jsd, intr;

    // will be modified
    Eigen::MatrixXd y;
    Eigen::VectorXd w;
    Eigen::MatrixXd a0;
    Eigen::VectorXd ao;
    Eigen::VectorXi ia;
    Eigen::VectorXi kin;
    Eigen::VectorXd rsqo; 
    Eigen::VectorXd almo;
    int nlp = 0, jerr = 0, lmu = 0;
    
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXi& jd;
    const Eigen::MatrixXd& cl;

    MultElnetDSCheckPack(
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
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha), flmin(_flmin), isd(_isd), jsd(_jsd), intr(_intr)
        , y(_y)
        , w(_w)
        , a0(_y.cols(), _nlam)
        , ao(_nx * _y.cols() * _nlam)
        , ia(nx)
        , kin(nlam)
        , rsqo(nlam)
        , almo(nlam)
        , X(_X)
        , ulam(_ulam), vp(_vp), jd(_jd)
        , cl(_cl)
    {
        ao.setZero();
        ia.setZero();
        kin.setZero();
        rsqo.setZero();
        almo.setZero();
    }

    void fit() 
    {
        Eigen::MatrixXd X_d = X;
        transl::multelnet<double>(
                alpha, X_d, y, w, jd, vp, cl, ne, nx, 
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, kin, rsqo, almo, nlp, jerr);
    }

    void fit_old() 
    {
        // apply my fix and check if the results match with dense
        transl::multspelnet<double>(
                alpha, X, y, w, jd, vp, cl, ne, nx, 
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, kin, rsqo, almo, nlp, jerr);
    }
};

struct multelnet_dscheck_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
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
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd y, cl;
    Eigen::VectorXd w, vp, ulam;
    Eigen::VectorXi jd;
    size_t maxit, nlam, nx, ne;
    double alpha, flmin;
    bool isd, jsd, intr;

    void check_pack(const MultElnetDSCheckPack& actual,
                    const MultElnetDSCheckPack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        expect_near_vec(actual.rsqo, expected.rsqo, 1e-14);

        expect_eq_vec(actual.kin, expected.kin);
        expect_eq_vec(actual.ia, expected.ia);

        EXPECT_EQ(actual.almo.size(), expected.almo.size());
        for (int i = 0; i < actual.almo.size(); ++i) {
            EXPECT_NEAR(actual.almo[i], expected.almo[i], actual.almo[i]*1e-15);
        }

        expect_near_mat(actual.a0, expected.a0, 1e-14);
        expect_near_vec(actual.ao, expected.ao, 1e-14);

        // Note: we do not compare y because dense and sparse treat it differently.
        // In dense, it is the full residual. 
        // In sparse, it's only a partial residual coming from the raw (uncentered) x_i
        // (i.e. it is missing the component that comes from centering of x_i).
    }
};

TEST_P(multelnet_dscheck_fixture, multelnet_dscheck_test)
{
    MultElnetDSCheckPack actual(
            maxit, nx, ne, nlam, isd, jsd, intr, alpha, flmin,
            X, y, w, ulam, vp, cl, jd);
    MultElnetDSCheckPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultElnetDSCheckSuite, multelnet_dscheck_fixture,
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
