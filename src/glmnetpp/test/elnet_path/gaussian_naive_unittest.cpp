#include <elnet_path/gaussian_base.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/elnet2.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_naive.hpp>

namespace glmnetpp {

struct GaussianNaivePack
    : GaussianPack
{
    Eigen::VectorXd y;
    const Eigen::MatrixXd& X;

    GaussianNaivePack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            const Eigen::MatrixXd& _X,
            Eigen::VectorXd& _y,
            Eigen::VectorXd& _xv,
            Eigen::VectorXd& _ulam,
            Eigen::VectorXd& _vp,
            Eigen::MatrixXd& _cl,
            Eigen::VectorXi& _ju)
        : GaussianPack(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _xv, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
    {}

    void fit() override
    {
        using internal_t = ElnetPointInternal<
                    util::glm_type::gaussian,
                    util::mode_type<util::glm_type::gaussian>::naive,
                    double,
                    int, 
                    int>;
        using elnet_point_t = ElnetPoint<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::naive,
                internal_t>;
        using elnet_path_t = ElnetPath<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::naive,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, y, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() override
    {
        transl::elnet2<double>(
                alpha, ju, vp, cl, y, ne, nx, X, nlam,
                flmin, ulam, thr, maxit, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }
};

struct gaussian_naive_fixture
    : gaussian_fixture
{
protected:
    using base_t = gaussian_fixture;
    void check_pack(const GaussianNaivePack& actual,
                    const GaussianNaivePack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_vec(actual.y, expected.y);
        expect_double_eq_mat(actual.X, expected.X);
    }
};

TEST_P(gaussian_naive_fixture, gaussian_naive_test)
{
    GaussianNaivePack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            X, y, xv, ulam, vp, cl, ju);
    GaussianNaivePack expected(actual);
    run(actual, expected, 4, 5); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    GaussianNaiveSuite, gaussian_naive_fixture,
    testing::Combine(
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
