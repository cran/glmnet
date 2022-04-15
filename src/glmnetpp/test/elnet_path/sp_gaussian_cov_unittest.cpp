#include <testutil/translation/spelnet1.hpp>
#include <testutil/mock_pb.hpp>
#include <elnet_path/gaussian_base.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_cov.hpp>

namespace glmnetpp {

struct SpGaussianCovPack
    : GaussianPack
{
    Eigen::VectorXd g;
    const Eigen::VectorXd& w;
    const Eigen::VectorXd& xm;
    const Eigen::VectorXd& xs;
    const Eigen::SparseMatrix<double>& X;

    SpGaussianCovPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            const Eigen::VectorXd& _w,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _xm,
            const Eigen::VectorXd& _xs,
            Eigen::VectorXd& _xv,
            Eigen::VectorXd& _ulam,
            Eigen::VectorXd& _vp,
            Eigen::MatrixXd& _cl,
            Eigen::VectorXi& _ju)
        : GaussianPack(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _xv, _ulam, _vp, _cl, _ju)
        , g(_X.transpose() * _y)
        , w(_w)
        , xm(_xm)
        , xs(_xs)
        , X(_X)
    {}

    void fit() override
    {
        using internal_t = SpElnetPointInternal<
                    util::glm_type::gaussian,
                    util::mode_type<util::glm_type::gaussian>::cov,
                    double,
                    int, 
                    int>;
        using elnet_point_t = SpElnetPoint<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::cov,
                internal_t>;
        using elnet_path_t = SpElnetPath<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::cov,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, g, w, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() override
    {
        transl::spelnet1<double>(
                alpha, ju, vp, cl, g, w, ne, nx, X, 
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, lmu, 
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }
};

struct sp_gaussian_cov_fixture
    : sp_gaussian_fixture
{
protected:
    using base_t = sp_gaussian_fixture;

    void check_pack(const SpGaussianCovPack& actual,
                    const SpGaussianCovPack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_vec(actual.g, expected.g);
    }
};

TEST_P(sp_gaussian_cov_fixture, sp_gaussian_cov_test)
{
    SpGaussianCovPack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            w, X, y, xm, xs, xv, ulam, vp, cl, ju);
    SpGaussianCovPack expected(actual);
    run(actual, expected, 6, 7);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpGaussianCovSuite, sp_gaussian_cov_fixture,
    testing::Combine(
        // seed, n, p, maxit, nlam, alpha, flmin
        testing::Values(241, 412, 23968, 31, 87201, 746, 104),
        testing::Values(10, 30, 50),
        testing::Values(5, 20, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 4, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5)
        )
);

} // namespace glmnetpp
