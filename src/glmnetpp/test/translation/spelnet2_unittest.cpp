#include <legacy/legacy.h>
#include <testutil/translation/spelnet2.hpp>
#include <translation/elnet_base_fixture.hpp>

namespace glmnetpp {

struct SpElnet2Pack
    : ElnetBasePack
{
    Eigen::VectorXd y;
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& w;
    const Eigen::VectorXd& xm;
    const Eigen::VectorXd& xs;

    SpElnet2Pack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _xm,
            const Eigen::VectorXd& _xs,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : ElnetBasePack(_maxit, _nx, _ne, _nlam, _alpha, _flmin, _xv, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , w(_w)
        , xm(_xm)
        , xs(_xs)
    {}

    void fit() override
    {
        transl::spelnet2<float>(
                alpha, ju, vp, cl, y, w, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();

        // for god-knows what reason, spelnet2 requires 1-indexed inner and outer indices :(
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);

        spelnet2_(const_cast<double*>(&alpha), &ni, y.data(),  
                  const_cast<double*>(w.data()), &no,
                  const_cast<int*>(&ne), 
                  const_cast<int*>(&nx),
                  const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(), 
                  const_cast<int*>(ju.data()), 
                  const_cast<double*>(vp.data()), 
                  const_cast<double*>(cl.data()), 
                  const_cast<int*>(&nlam),
                  const_cast<double*>(&flmin), 
                  const_cast<double*>(ulam.data()), 
                  const_cast<double*>(&thr), 
                  const_cast<int*>(&maxit), 
                  const_cast<double*>(xm.data()), 
                  const_cast<double*>(xs.data()),
                  const_cast<double*>(xv.data()), &lmu,
                  ao.data(), ia.data(), kin.data(), rsqo.data(), almo.data(),
                  &nlp, &jerr);
    }
};

struct spelnet2_fixture
    : spelnet_base_fixture
{
protected:
    using base_t = spelnet_base_fixture;

    void check_pack(const SpElnet2Pack& actual,
                    const SpElnet2Pack& expected)
    {
        base_t::check_pack(actual, expected);
        expect_float_eq_vec(actual.y, expected.y);
    }
};

TEST_P(spelnet2_fixture, spelnet2_test)
{
    SpElnet2Pack actual(
            maxit, nx, ne, nlam, alpha, flmin,
            w, xm ,xs, X, y, xv, ulam, vp, cl, ju);
    SpElnet2Pack expected(actual);
    run(actual, expected, 10, 11); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpElnet2Suite, spelnet2_fixture,
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
