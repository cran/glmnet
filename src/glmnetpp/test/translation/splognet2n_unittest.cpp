#include <legacy/legacy.h>
#include <testutil/translation/splognet2n.hpp>
#include <translation/lognet_base_fixture.hpp>

namespace glmnetpp {

struct SpLognet2nPack: LognetBasePack
{
    const Eigen::VectorXd& y;
    const Eigen::SparseMatrix<double>& X;
    Eigen::VectorXd g;  // g is offset
    const Eigen::VectorXd& xb, xs;

    SpLognet2nPack(
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
        : LognetBasePack(_maxit, _nx, _ne, _nlam, _kopt, _isd, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
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
        transl::splognet2n<float, false>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, xb, xs, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        sprlognet2n_(
                const_cast<double*>(&alpha), &no, &ni, 
                const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(), 
                const_cast<double*>(y.data()),
                g.data(),
                const_cast<double*>(w.data()),
                const_cast<int*>(ju.data()), 
                const_cast<double*>(vp.data()), 
                const_cast<double*>(cl.data()),
                const_cast<int*>(&ne), 
                const_cast<int*>(&nx), 
                const_cast<int*>(&nlam),
                const_cast<double*>(&flmin), 
                const_cast<double*>(ulam.data()), 
                const_cast<double*>(&thr), 
                const_cast<int*>(&isd),
                const_cast<int*>(&intr),
                const_cast<int*>(&maxit), 
                const_cast<int*>(&kopt),
                const_cast<double*>(xb.data()), 
                const_cast<double*>(xs.data()),
                &lmu,
                a0.data(), a.data(), m.data(), kin.data(), &dev0, dev.data(), alm.data(),
                &nlp, &jerr);
    }
};

struct splognet2n_fixture
    : splognet_base_fixture
{
protected:
    using base_t = lognet_base_fixture;

    void check_pack(const SpLognet2nPack& actual,
                    const SpLognet2nPack& expected)
    {
        base_t::check_pack(actual, expected);    
        expect_near_vec(actual.g, expected.g, 1e-7);
    }
};

TEST_P(splognet2n_fixture, splognet2n_test)
{
    SpLognet2nPack actual(
            maxit, nx, ne, nlam, alpha, flmin, kopt, isd, intr,
            X, y, w, ulam, vp, cl, ju, xm, xs);
    SpLognet2nPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpLognet2nSuite, splognet2n_fixture,
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
