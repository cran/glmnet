#include <legacy/legacy.h>
#include <testutil/translation/lognet2n.hpp>
#include <translation/lognet_base_fixture.hpp>

namespace glmnetpp {

struct Lognet2nPack: LognetBasePack
{
    const Eigen::VectorXd& y;
    const Eigen::MatrixXd& X;
    Eigen::VectorXd g;  // g is offset

    Lognet2nPack(
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
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : LognetBasePack(_maxit, _nx, _ne, _nlam, _kopt, _isd, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , g(_g)
    {}

    void fit() override
    {
        transl::lognet2n<float>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        lognet2n_(const_cast<double*>(&alpha), &no, &ni, 
                const_cast<double*>(X.data()), 
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
                &lmu,
                a0.data(), a.data(), m.data(), kin.data(), &dev0, dev.data(), alm.data(),
                &nlp, &jerr);
    }
};

struct lognet2n_fixture
    : lognet_base_fixture
{
protected:
    using base_t = lognet_base_fixture;

    void check_pack(const Lognet2nPack& actual,
                    const Lognet2nPack& expected)
    {
        base_t::check_pack(actual, expected);    
        expect_near_vec(actual.g, expected.g, 1e-7);
    }
};

TEST_P(lognet2n_fixture, lognet2n_test)
{
    Lognet2nPack actual(
            maxit, nx, ne, nlam, alpha, flmin, kopt, isd, intr,
            X, y, g, w, ulam, vp, cl, ju);
    Lognet2nPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    Lognet2nSuite, lognet2n_fixture,
    testing::Combine(
        testing::Values(3029,238,482873,3887472,1039284,38888328), // seed
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
