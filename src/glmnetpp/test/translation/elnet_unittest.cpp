#include <legacy/legacy.h>
#include <testutil/translation/elnet.hpp>
#include <translation/elnet_base_fixture.hpp>

namespace glmnetpp {

struct ElnetPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, intr, ka;

    Eigen::MatrixXd X;
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

    ElnetPack(
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
        transl::elnet<float>(
                ka, alpha, X, y, w, jd, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                a0, ca, ia, nin, rsq, alm, nlp, jerr);
    }

    void fit_old()
    {
        int no = X.rows();
        int ni = X.cols();
        int iisd = isd;
        int iintr = intr;
        int ika = ka + 1;
        ::elnet_(
                &ika, 
                const_cast<double*>(&alpha), &no, &ni,
                X.data(), y.data(), w.data(),
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
                nin.data(), rsq.data(), alm.data(), &nlp, &jerr);
    }
};

struct elnet_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X(n, p);
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
    Eigen::MatrixXd X, cl;
    Eigen::VectorXd y, w, ulam, vp;
    Eigen::VectorXi jd;
    int nx, ne, maxit, nlam;
    double alpha, flmin;
    bool isd, intr, ka;

    void check_pack(const ElnetPack& actual,
                    const ElnetPack& expected)
    {
        expect_float_eq_mat(actual.X, expected.X);
        expect_float_eq_vec(actual.w, expected.w);

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

        expect_near_vec(actual.ia, expected.ia, 1);
        expect_eq_vec(actual.nin, expected.nin);
        expect_float_eq_mat(actual.ca, expected.ca);
        expect_float_eq_vec(actual.a0, expected.a0);
        expect_float_eq_vec(actual.y, expected.y);
        expect_float_eq_vec(actual.rsq, expected.rsq);
        expect_float_eq_vec(actual.alm, expected.alm);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
    }
};

TEST_P(elnet_fixture, elnet_test)
{
    ElnetPack actual(
            maxit, nx, ne, nlam, alpha, flmin, isd, intr, ka,
            X, y, w, ulam, vp, cl, jd);
    ElnetPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    ElnetSuite, elnet_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin, isd, intr, ka
        testing::Values(241, 412, 23968),
        testing::Values(10, 30, 50),
        testing::Values(5, 20, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 4, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5),
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

} // namespace glmnetpp
