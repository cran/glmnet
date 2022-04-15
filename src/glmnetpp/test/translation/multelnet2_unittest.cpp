#include <legacy/legacy.h>
#include <testutil/translation/multelnet2.hpp>
#include <translation/elnet_base_fixture.hpp>

namespace glmnetpp {

struct MultElnet2Pack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;

    // will be modified
    Eigen::MatrixXd y;
    Eigen::VectorXd ao;
    Eigen::VectorXi ia;
    Eigen::VectorXi kin;
    Eigen::VectorXd rsqo; 
    Eigen::VectorXd almo;
    double ys0;
    int nlp = 0, jerr = 0, lmu = 0;
    
    const Eigen::MatrixXd& X;
    const Eigen::VectorXd& xv;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXd& cl;
    const Eigen::VectorXi& ju;

    MultElnet2Pack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            double _ys0,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::VectorXd& _cl,
            const Eigen::VectorXi& _ju)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha), flmin(_flmin)
        , y(_y)
        , ao(_nx * _y.cols() * _nlam)
        , ia(nx)
        , kin(nlam)
        , rsqo(nlam)
        , almo(nlam)
        , ys0(_ys0)
        , X(_X), xv(_xv), ulam(_ulam), vp(_vp), cl(_cl), ju(_ju)
    {
        ao.setZero();
        ia.setZero();
        kin.setZero();
        rsqo.setZero();
        almo.setZero();
    }

    void fit()
    {
        transl::multelnet2<float>(
                alpha, ju, vp, cl, y, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xv, ys0, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }

    void fit_old() 
    {
        int ni = X.cols();
        int nr = y.cols();
        int no = X.rows();
        multelnet2_(
                const_cast<double*>(&alpha), &ni, &nr,
                const_cast<int*>(ju.data()), 
                const_cast<double*>(vp.data()), 
                const_cast<double*>(cl.data()),
                const_cast<double*>(y.data()), &no, 
                const_cast<int*>(&ne), 
                const_cast<int*>(&nx), 
                const_cast<double*>(X.data()), 
                const_cast<int*>(&nlam),
                const_cast<double*>(&flmin), 
                const_cast<double*>(ulam.data()), 
                const_cast<double*>(&thr), 
                const_cast<int*>(&maxit), 
                const_cast<double*>(xv.data()), &ys0, &lmu,
                ao.data(), ia.data(), kin.data(), rsqo.data(), almo.data(),
                &nlp, &jerr);
    }
};

struct multelnet2_fixture
    : elnet_base_fixture
{
    void SetUp() override
    {
        size_t seed, n, p, nr = 3;
        bool isd = true, jsd = true, intr = true;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p, 0.85);
        y.resize(n, nr);
        for (size_t i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        auto w = dgen.make_w(n);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl.resize(2 * nr * p);
        Eigen::MatrixXd cl_sub = dgen.make_cl(p);
        // TODO
        cl_sub.row(0).array() = -0.5;
        cl_sub.row(1).array() = 0.5;
        for (size_t i = 0; i < p; ++i) {
            Eigen::Map<Eigen::MatrixXd> cl_slice(
                    cl.data() + i * 2 * nr, 2, nr);
            cl_slice.colwise() = cl_sub.col(i);
        }
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        xv.setOnes(p);

        Eigen::VectorXd xm(p), xs(p), ym(nr), ys(nr);
        MultStandardize1::eval(X, y, w, isd, jsd, intr, ju, xm, xs, ym, ys, xv, ys0);
    }

protected:
    using base_t = elnet_base_fixture;
    Eigen::MatrixXd y;
    Eigen::VectorXd cl;
    double ys0;

    void check_pack(const MultElnet2Pack& actual,
                    const MultElnet2Pack& expected)
    {
        expect_near_vec(actual.ao, expected.ao, 1e-7);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.kin, expected.kin);
        expect_float_eq_vec(actual.rsqo, expected.rsqo);
        expect_float_eq_vec(actual.almo, expected.almo);
        EXPECT_FLOAT_EQ(actual.ys0, expected.ys0);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        if (actual.jerr == util::bnorm_maxit_reached_error().err_code())
            std::cout << actual.jerr << std::endl;
        EXPECT_EQ(actual.lmu, expected.lmu);
        if (actual.jerr > util::max_active_reached_error().err_code(0)) {
            expect_near_mat(actual.y, expected.y, 1e-7);
        }
    }
};

TEST_P(multelnet2_fixture, multelnet2_test)
{
    MultElnet2Pack actual(
            maxit, nx, ne, nlam, alpha, flmin, ys0,
            X, y, xv, ulam, vp, cl, ju);
    MultElnet2Pack expected(actual);
    run(actual, expected, 4, 5);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultElnet2Suite, multelnet2_fixture,
    testing::Combine(
        testing::Values(3827, 1039298428, 2, 314, 2985), // seed
        testing::Values(10, 30, 50, 1000),          // n
        testing::Values(5, 20, 40, 100),           // p
        testing::Values(1, 50, 100000),          // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5)       // flmin
        )
);

} // namespace glmnetpp
