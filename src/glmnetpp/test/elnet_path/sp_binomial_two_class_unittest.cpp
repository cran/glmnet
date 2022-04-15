#include <elnet_path/binomial_base.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/splognet2n.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp> // needs definition of naive base
#include <glmnetpp_bits/elnet_point/binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_point/sp_binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_point/sp_binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/sp_binomial_two_class.hpp>

namespace glmnetpp {

struct SpBinomialTwoClassPack: BinomialPack
{
    const Eigen::VectorXd& y;
    const Eigen::SparseMatrix<double>& X;
    Eigen::VectorXd g;  // g is offset
    const Eigen::VectorXd& xb, xs;

    SpBinomialTwoClassPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            int _kopt,
            int _isd,
            int _intr,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju,
            const Eigen::VectorXd& _xb,
            const Eigen::VectorXd& _xs)
        : BinomialPack(_maxit, _nx, _ne, _nlam, _kopt, _isd, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , g(_y.size())
        , xb(_xb)
        , xs(_xs)
    {
        g.setOnes(); 
    }

    void fit() override
    {
        using internal_t = SpElnetPointInternal<
                    util::glm_type::binomial,
                    util::mode_type<util::glm_type::binomial>::two_class,
                    double,
                    int, 
                    int>;
        using elnet_point_t = SpElnetPoint<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::two_class,
                internal_t>;
        using elnet_path_t = SpElnetPath<
            util::glm_type::binomial,
            util::mode_type<util::glm_type::binomial>::two_class,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, ne, nx, X, y, g, w, 
                nlam, flmin, ulam, xb, xs, thr, isd, intr, maxit, kopt, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() override
    {
        transl::splognet2n<double, true>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, xb, xs, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }
};

struct sp_binomial_two_class_fixture
    : sp_binomial_fixture
{
protected:
    using base_t = sp_binomial_fixture;

    void check_pack(const SpBinomialTwoClassPack& actual,
                    const SpBinomialTwoClassPack& expected)
    {
        base_t::check_pack(actual, expected);    
        expect_near_vec(actual.g, expected.g, 1e-7);
    }
};

TEST_P(sp_binomial_two_class_fixture, sp_binomial_two_class_test)
{
    SpBinomialTwoClassPack actual(
            maxit, nx, ne, nlam, alpha, flmin, kopt, isd, intr,
            X, y, w, ulam, vp, cl, ju, xm, xs);
    SpBinomialTwoClassPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpBinomialTwoClassSuite, sp_binomial_two_class_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5),       // flmin
        testing::Bool(),                      // isd
        testing::Bool(),                      // intr
        testing::Values(0,1,2)                // kopt
        )
);

} // namespace glmnetpp
