#include <testutil/translation/elnet1.hpp>
#include <testutil/mock_pb.hpp>
#include <elnet_path/gaussian_base.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_cov.hpp>

namespace glmnetpp {

struct GaussianCovPack
    : GaussianPack
{
    Eigen::VectorXd g;
    const Eigen::MatrixXd& X;

    GaussianCovPack(int _maxit,
                    int _nx,
                    int _ne,
                    int _nlam,
                    double _alpha,
                    double _flmin,
                    const Eigen::MatrixXd& _X,
                    const Eigen::VectorXd& _y,
                    Eigen::VectorXd& _xv,
                    Eigen::VectorXd& _ulam,
                    Eigen::VectorXd& _vp,
                    Eigen::MatrixXd& _cl,
                    Eigen::VectorXi& _ju)
        : GaussianPack(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _xv, _ulam, _vp, _cl, _ju)
        , g(_X.transpose() * _y)
        , X(_X)
    {}

    void fit() override
    {
        using internal_t = ElnetPointInternal<
                    util::glm_type::gaussian,
                    util::mode_type<util::glm_type::gaussian>::cov,
                    double,
                    int, 
                    int>;
        using elnet_point_t = ElnetPoint<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::cov,
                internal_t>;
        using elnet_path_t = ElnetPath<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::cov,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, g, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() override
    {
        transl::elnet1<double>(
                alpha, ju, vp, cl, g, ne, nx, X, 
                nlam, flmin, ulam, thr, maxit, xv, lmu, 
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }
};

struct gaussian_cov_fixture
    : gaussian_fixture
{
protected:
    using base_t = gaussian_fixture;

    void check_pack(const GaussianCovPack& actual,
                    const GaussianCovPack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_near_vec(actual.g, expected.g, 1e-15);
        expect_double_eq_mat(actual.X, expected.X);
    }
};

TEST_P(gaussian_cov_fixture, gaussian_cov_test)
{
    GaussianCovPack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            X, y, xv, ulam, vp, cl, ju);
    GaussianCovPack expected(actual);
    run(actual, expected, 0, 1); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    GaussianCovSuite, gaussian_cov_fixture,
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
