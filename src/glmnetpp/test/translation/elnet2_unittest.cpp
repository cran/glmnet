#include <legacy/legacy.h>
#include <testutil/translation/elnet2.hpp>
#include <translation/elnet_base_fixture.hpp>

namespace glmnetpp {

struct Elnet2Pack: ElnetBasePack
{
    Eigen::VectorXd y;
    const Eigen::MatrixXd& X;

    Elnet2Pack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            Eigen::MatrixXd& _X,
            const Eigen::VectorXd& _y,
            Eigen::VectorXd& _xv,
            Eigen::VectorXd& _ulam,
            Eigen::VectorXd& _vp,
            Eigen::MatrixXd& _cl,
            Eigen::VectorXi& _ju)
        : ElnetBasePack(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _xv, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
    {}

    void fit() override
    {
        transl::elnet2<float>(
                alpha, ju, vp, cl, y, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        elnet2_(const_cast<double*>(&alpha), &ni, 
                const_cast<int*>(ju.data()), 
                const_cast<double*>(vp.data()), 
                const_cast<double*>(cl.data()),
                y.data(), &no, 
                const_cast<int*>(&ne), 
                const_cast<int*>(&nx), 
                const_cast<double*>(X.data()), 
                const_cast<int*>(&nlam),
                const_cast<double*>(&flmin), 
                const_cast<double*>(ulam.data()), 
                const_cast<double*>(&thr), 
                const_cast<int*>(&maxit), 
                const_cast<double*>(xv.data()), &lmu,
                ao.data(), ia.data(), kin.data(), rsqo.data(), almo.data(),
                &nlp, &jerr);
    }
};

struct elnet2_fixture
    : elnet_base_fixture
{
protected:
    using base_t = elnet_base_fixture;

    void check_pack(const Elnet2Pack& actual,
                    const Elnet2Pack& expected)
    {
        base_t::check_pack(actual, expected);    
        expect_float_eq_vec(actual.y, expected.y);
    }
};

TEST_P(elnet2_fixture, elnet2_test)
{
    Elnet2Pack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            X, y, xv, ulam, vp, cl, ju);
    Elnet2Pack expected(actual);
    run(actual, expected, 4, 5);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    Elnet2Suite, elnet2_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5)        // flmin
        )
);

} // namespace glmnetpp
