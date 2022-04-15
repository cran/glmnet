#include <legacy/legacy.h>
#include <testutil/translation/lognetn.hpp>
#include <translation/lognet_base_fixture.hpp>

namespace glmnetpp {

struct LognetnPack: LognetBasePack
{
    const Eigen::MatrixXd& y;
    const Eigen::MatrixXd& X;
    const int nc;
    Eigen::MatrixXd g;  // g is offset
    Eigen::MatrixXd a0;
    Eigen::VectorXd a;  // flattened 3-d array

    LognetnPack(
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
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _ju)
        : LognetBasePack(_maxit, _nx, _ne, _nlam, _kopt, _isd, _intr, _alpha, _flmin, _w, _ulam, _vp, _cl, _ju)
        , y(_y)
        , X(_X)
        , nc(_y.cols())
        , g(_y.rows(), _y.cols())
        , a0(_y.cols(), _nlam)
        , a(_nx * _y.cols() * _nlam)
    {
        g.setOnes(); 
    }

    void fit() override
    {
        transl::lognetn<float>(
                alpha, X, y, g, w, ju, vp, cl, ne, nx,
                nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                a0, a, m, kin, dev0, dev, alm, nlp, jerr);
    }

    void fit_old() override
    {
        int ni = X.cols();
        int no = X.rows();
        int nc = y.cols();
        lognetn_(const_cast<double*>(&alpha), &no, &ni, &nc,
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

struct lognetn_fixture
    : lognet_base_fixture
{
    void SetUp() override
    {
        base_t::SetUp();
        y.resize(base_t::y.rows(), 2);
        y.col(0) = base_t::y;
        y.col(1).array() = 1-base_t::y.array();
    }

protected:
    using base_t = lognet_base_fixture;
    Eigen::MatrixXd y;

    void check_pack(const LognetnPack& actual,
                    const LognetnPack& expected)
    {
        base_t::check_pack(actual, expected);    
        expect_near_mat(actual.g, expected.g, 1e-7);

        bool max_active_reached = (actual.jerr <= util::max_active_reached_error().err_code(0));

        // Note: the last column may not be kosher to compare with because fortran has a bug.
        int cols = actual.a0.cols() - max_active_reached;
        if (cols >= 0) {
            expect_near_mat(actual.a0.block(actual.a0.rows(), cols, 0, 0), 
                            expected.a0.block(expected.a0.rows(), cols, 0, 0), 1e-7);
        }
        int nx = actual.nx;
        int nc = actual.nc;
        int nlam = actual.nlam - max_active_reached;
        // Same thing: last slice may not be kosher to compare with because fortran has a bug.

        for (int i = 0; i < nlam; ++i) {
            Eigen::Map<const Eigen::MatrixXd> actual_a_slice(
                    actual.a.data() + nx * nc * i, nx, nc);
            Eigen::Map<const Eigen::MatrixXd> expected_a_slice(
                    expected.a.data() + nx * nc * i, nx, nc);
            expect_near_mat(actual_a_slice, expected_a_slice, 1e-7);
        }
    }
};

TEST_P(lognetn_fixture, lognetn_test)
{
    LognetnPack actual(
            maxit, nx, ne, nlam, alpha, flmin, kopt, isd, intr,
            X, y, w, ulam, vp, cl, ju);
    LognetnPack expected(actual);
    run(actual, expected, 0, 1);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    LognetnSuite, lognetn_fixture,
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
