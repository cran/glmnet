#include <elnet_path/binomial_base.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/multsplognetn.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_binomial_multi_class_group.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp> // for naive and multi base
#include <glmnetpp_bits/elnet_point/sp_binomial_multi_class_group.hpp>
#include <glmnetpp_bits/elnet_point/binomial_multi_class_group.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/sp_binomial_multi_class_group.hpp>

namespace glmnetpp {

struct SpBinomialMultiClassGroupPack: BinomialPack
{
    const Eigen::MatrixXd& y;
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& xv;
    const Eigen::VectorXd& xm;
    const Eigen::VectorXd& xs;
    const int nc;
    Eigen::MatrixXd g;  // g is offset
    Eigen::MatrixXd a0;
    Eigen::VectorXd a;  // flattened 3-d array

    SpBinomialMultiClassGroupPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            int _intr,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _xm,
            const Eigen::VectorXd& _xs,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : BinomialPack(_maxit, _nx, _ne, _nlam, 2, true, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , xv(_xv)
        , xm(_xm)
        , xs(_xs)
        , nc(_y.cols())
        , g(_y.rows(), _y.cols())
        , a0(_y.cols(), _nlam)
        , a(_nx * _y.cols() * _nlam)
    {
        g.setZero(); 
        a0.setZero();
        a.setZero();
    }

    void fit() override
    {
        using internal_t = SpElnetPointInternal<
                    util::glm_type::binomial,
                    util::mode_type<util::glm_type::binomial>::multi_class_group,
                    double,
                    int, 
                    int>;
        using elnet_point_t = SpElnetPoint<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::multi_class_group,
                internal_t>;
        using elnet_path_t = SpElnetPath<
            util::glm_type::binomial,
            util::mode_type<util::glm_type::binomial>::multi_class_group,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, ne, nx, X, y, g, w, 
                nlam, flmin, ulam, thr, intr, maxit, xm, xs, xv, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() override
    {
        transl::multsplognetn<double, true>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, intr, maxit, xv, xm, xs, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }
};

struct sp_binomial_multi_class_group_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        int isd;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X_sparse(n, p);
        auto beta = dgen.make_beta(p);
        y.resize(n, 2);
        y.col(0) = dgen.make_y(X, beta);
        y.col(0).array() =  (y.col(0).array() > y.col(0).mean()).template cast<double>();
        y.col(1).array() = 1.-y.col(0).array();
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);

        xm.setZero(p);
        xs.setZero(p);
        xv.setZero(p);
        MultSpLStandardize2::eval(X, w, ju, isd, intr, xm, xs, xv);
    }

protected:
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd cl, y;
    Eigen::VectorXd xv, ulam, vp, w, xm, xs;
    Eigen::VectorXi ju;
    int nx, ne, maxit, nlam, intr;
    double alpha, flmin;

    void check_pack(const SpBinomialMultiClassGroupPack& actual,
                    const SpBinomialMultiClassGroupPack& expected)
    {
        expect_near_mat(actual.a0, expected.a0, 1e-15);
        expect_near_vec(actual.a, expected.a, 1e-15);
        if (actual.jerr > util::max_active_reached_error().err_code(0)) {
            expect_near_mat(actual.g, expected.g, 2e-15);
        }
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_DOUBLE_EQ(actual.dev0, expected.dev0);
        expect_near_vec(actual.dev, expected.dev, 1e-15);
        expect_eq_vec(actual.kin, expected.kin);
        expect_eq_vec(actual.m, expected.m);

        EXPECT_EQ(actual.alm.size(), expected.alm.size());
        for (int i = 0; i < actual.alm.size(); ++i) {
            EXPECT_NEAR(actual.alm[i], expected.alm[i], actual.alm[i]*1e-15);
        }
    }
};

TEST_P(sp_binomial_multi_class_group_fixture, sp_binomial_multi_class_group_test)
{
    SpBinomialMultiClassGroupPack actual(
            maxit, nx, ne, nlam, alpha, flmin, intr,
            X, y, w, ulam, xv, xm, xs, vp, cl, ju);
    SpBinomialMultiClassGroupPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpBinomialMultiClassGroupSuite, sp_binomial_multi_class_group_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31, 1242, 4552), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5),       // flmin
        testing::Bool(),                      // isd
        testing::Bool()                       // intr
        )
);

} // namespace glmnetpp
