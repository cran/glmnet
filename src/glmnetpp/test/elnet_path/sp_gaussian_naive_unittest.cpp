#include <testutil/translation/spelnet2.hpp>
#include <testutil/mock_pb.hpp>
#include <elnet_path/gaussian_base.hpp>
#include <glmnetpp_bits/internal.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_naive.hpp>

namespace glmnetpp {

struct SpGaussianNaivePack
    : GaussianPack
{
    Eigen::VectorXd y;
    const Eigen::VectorXd& w;
    const Eigen::VectorXd& xm;
    const Eigen::VectorXd& xs;
    const Eigen::SparseMatrix<double>& X;

    SpGaussianNaivePack(
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
        , y(_y)
        , w(_w)
        , xm(_xm)
        , xs(_xs)
        , X(_X)
    {}

    void fit() override
    {
        using internal_t = SpElnetPointInternal<
                    util::glm_type::gaussian,
                    util::mode_type<util::glm_type::gaussian>::naive,
                    double,
                    int, 
                    int>;
        using elnet_point_t = SpElnetPoint<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::naive,
                internal_t>;
        using elnet_path_t = SpElnetPath<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::naive,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, y, w, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_transl() override
    {
        transl::spelnet2<double>(
                alpha, ju, vp, cl, y, w, ne, nx, X, 
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, lmu, 
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }
};

struct sp_gaussian_naive_fixture
    : sp_gaussian_fixture
{
protected:
    using base_t = sp_gaussian_fixture;

    void check_pack(const SpGaussianNaivePack& actual,
                    const SpGaussianNaivePack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_double_eq_vec(actual.y, expected.y);
    }
};

TEST_P(sp_gaussian_naive_fixture, sp_gaussian_naive_test)
{
    SpGaussianNaivePack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            w, X, y, xm, xs, xv, ulam, vp, cl, ju);
    SpGaussianNaivePack expected(actual);
    run(actual, expected, 12, 13);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpGaussianNaiveSuite, sp_gaussian_naive_fixture,
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
