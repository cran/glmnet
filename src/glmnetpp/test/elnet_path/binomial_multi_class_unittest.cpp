#include <elnet_path/binomial_base.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/lognetn.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/internal/binomial_multi_class.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp> // naive
#include <glmnetpp_bits/elnet_point/binomial_multi_class.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/binomial_multi_class.hpp>

namespace glmnetpp {

struct BinomialMultiClassPack: BinomialPack
{
    const Eigen::MatrixXd& y;
    const Eigen::MatrixXd& X;
    const int nc;
    Eigen::MatrixXd g;  // g is offset

    // override base class members
    Eigen::MatrixXd a0;
    Eigen::VectorXd a;  // flattened 3-d array

    BinomialMultiClassPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            int _kopt,
            int _isd,
            int _intr,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : BinomialPack(_maxit, _nx, _ne, _nlam, _kopt, _isd, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , nc(_y.cols())
        , g(_g)
        , a0(_y.cols(), _nlam)
        , a(_nx * _y.cols() * _nlam)
    {
        a0.setZero();
        a.setZero();
    }

    void fit() override
    {
        using internal_t = ElnetPointInternal<
                    util::glm_type::binomial,
                    util::mode_type<util::glm_type::binomial>::multi_class,
                    double,
                    int, 
                    int>;
        using elnet_point_t = ElnetPoint<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::multi_class,
                internal_t>;
        using elnet_path_t = ElnetPath<
            util::glm_type::binomial,
            util::mode_type<util::glm_type::binomial>::multi_class,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, ne, nx, X, y, g, w, 
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() override
    {
        transl::lognetn<double>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }
};

struct binomial_multi_class_fixture
    : binomial_fixture
{
    void SetUp() override
    {
        base_t::SetUp();
        y.resize(base_t::y.rows(), 2);
        y.col(0) = base_t::y;
        y.col(1).array() = 1-base_t::y.array();
    }

protected:
    using base_t = binomial_fixture;
    Eigen::MatrixXd y;

    void check_pack(const BinomialMultiClassPack& actual,
                    const BinomialMultiClassPack& expected)
    {
        base_t::check_pack(actual, expected);    
        expect_near_mat(actual.a0, expected.a0, 1e-15);
        expect_near_vec(actual.a, expected.a, 1e-15);
        expect_near_mat(actual.g, expected.g, 2e-15);
    }
};

TEST_P(binomial_multi_class_fixture, binomial_multi_class_test)
{
    BinomialMultiClassPack actual(
            maxit, nx, ne, nlam, alpha, flmin, kopt, isd, intr,
            X, y, w, g, ulam, vp, cl, ju);
    BinomialMultiClassPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    BinomialMultiClassSuite, binomial_multi_class_fixture,
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
