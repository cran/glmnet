#include <legacy/legacy.h>
#include <testutil/translation/elnet1.hpp>
#include <translation/elnet_base_fixture.hpp>

namespace glmnetpp {

struct Elnet1Pack
    : ElnetBasePack
{
    Eigen::VectorXd g;
    const Eigen::MatrixXd& X;

    Elnet1Pack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            const Eigen::MatrixXd& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : ElnetBasePack(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _xv, _ulam, _vp, _cl, _ju)
        , g(_X.transpose() * _y)
        , X(_X)
    {}

    void fit() override
    {
        transl::elnet1<float>(
                alpha, ju, vp, cl, g, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        elnet1_(const_cast<double*>(&alpha), &ni, 
                const_cast<int*>(ju.data()), 
                const_cast<double*>(vp.data()), 
                const_cast<double*>(cl.data()),
                g.data(), &no, 
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

struct elnet1_fixture
    : elnet_base_fixture
{
protected:
    using base_t = elnet_base_fixture;

    void check_pack(const Elnet1Pack& actual,
                    const Elnet1Pack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_float_eq_vec(actual.g, expected.g);
    }
};

TEST_P(elnet1_fixture, elnet1_test)
{
    Elnet1Pack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            X, y, xv, ulam, vp, cl, ju);
    Elnet1Pack expected(actual);
    run(actual, expected, 2, 3);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    Elnet1Suite, elnet1_fixture,
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
