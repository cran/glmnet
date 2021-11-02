#include <testutil/base_fixture.hpp>
#include <glmnetpp_bits/wls.hpp>
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <testutil/thread.hpp>

namespace glmnetpp {

struct WLSPack 
{
    // these values either really won't impact the algorithm
    // or just needs to be default-initialized
    int jerr = 0;
    int nlp = 0;
    int nino = 0;
    int m = 1242;
    double rsqc = 0;
    double thr = 1e-14;
    int iz = 1;

    // initialized upon construction
    double aint;
    int nx;
    const Eigen::MatrixXd& X, cl;
    const Eigen::VectorXd& vp, w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd beta, y, a, g, r, xv;
    Eigen::VectorXi ia ,iy, mm;
    double alm0;
    double almc;
    double alpha;
    int intr;
    int maxit;

    WLSPack(double _aint,
            int _nx,
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXd& _vp,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            const Eigen::VectorXd& _beta,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _a,
            const Eigen::VectorXd& _g,
            const Eigen::VectorXd& _r,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXi& _ia,
            const Eigen::VectorXi& _iy,
            const Eigen::VectorXi& _mm,
            double _alm0,
            double _almc,
            double _alpha,
            int _intr,
            int _maxit)
        : aint(_aint), nx(_nx), X(_X), cl(_cl), vp(_vp), w(_w), ju(_ju)
        , beta(_beta), y(_y), a(_a), g(_g), r(_r), xv(_xv)
        , ia(_ia), iy(_iy), mm(_mm)
        , alm0(_alm0)
        , almc(_almc)
        , alpha(_alpha)
        , intr(_intr)
        , maxit(_maxit)
    {}

    void run_wls() {
        int n = X.rows();
        int p = X.cols();
        wls(alm0, almc, alpha, m, n, p, X, r, xv, 
            w, intr, ju, vp, 
            cl, nx, thr, maxit, a, aint, g,
            ia, iy, iz, mm, nino, rsqc, nlp, jerr);
    }

    void run_wls_() {
        int n = X.rows();
        int p = X.cols();
        wls_(&alm0, &almc, &alpha, &m, &n, &p, 
             const_cast<double*>(X.data()), r.data(), 
             const_cast<double*>(w.data()), &intr, 
             const_cast<int*>(ju.data()), 
             const_cast<double*>(vp.data()), 
             const_cast<double*>(cl.data()), 
             &nx, &thr, &maxit, a.data(), &aint, g.data(),
             ia.data(), iy.data(), &iz, mm.data(), &nino, &rsqc, &nlp, &jerr);
    }
};

struct wls_fixture : 
    base_fixture,
    ::testing::WithParamInterface<
        std::tuple<int, int, int, int, int>
    >
{
protected:

    void test_wls_pack(const WLSPack& actual,
                       const WLSPack& expected)
    {
        expect_double_eq_vec(actual.r, expected.r);
        expect_double_eq_vec(actual.a, expected.a);
        EXPECT_DOUBLE_EQ(actual.aint, expected.aint);
        expect_double_eq_vec(actual.g, expected.g);
        expect_double_eq_vec(actual.ia, expected.ia);
        expect_double_eq_vec(actual.iy, expected.iy);
        EXPECT_EQ(actual.iz, expected.iz);
        expect_double_eq_vec(actual.mm, expected.mm);
        EXPECT_EQ(actual.nino, expected.nino);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_DOUBLE_EQ(actual.rsqc, expected.rsqc);
    }
        
};

TEST_P(wls_fixture, wls_test)
{
    double alm0 = 0.002;
    double almc = 0.001;
    double alpha = 0.9;

    int seed, n, p;
    bool intr, maxit;
    std::tie(seed, n, p, intr, maxit) = GetParam();
    DataGen dgen(seed);
    auto X = dgen.make_X(n, p);
    auto cl = dgen.make_cl(p);
    auto vp = dgen.make_vp(p);
    auto weights = dgen.make_w(n);
    auto inclusion = dgen.make_ju(p);
    auto beta = dgen.make_beta(p);
    auto y = dgen.make_y(X, beta);
    double aint = 0;
    auto nx = dgen.make_nx(p);

    Eigen::VectorXd a(p); a.setZero();
    Eigen::VectorXd g = X.transpose() * y;
    Eigen::VectorXd r(n); r.setZero();
    Eigen::VectorXd xv(p); xv.setOnes();
    Eigen::VectorXi ia(p); ia.setZero();
    Eigen::VectorXi iy(p); iy.setZero();
    Eigen::VectorXi mm(p); mm.setZero();

    for (int i = 0; i < std::min(p, 10); ++i) {
        ia(i) = mm(i) = i+1;
    }

    WLSPack actual_pack(
            aint, nx, X, cl, vp, weights, inclusion, beta, y,
            a, g, r, xv, ia, iy, mm,
            alm0, almc, alpha, intr, maxit);
    WLSPack expected_pack(actual_pack);

    std::thread actual_thr([&]() { actual_pack.run_wls(); });
    std::thread expected_thr([&]() { expected_pack.run_wls_(); });

    set_affinity(0, actual_thr.native_handle());
    set_affinity(1, expected_thr.native_handle());

    actual_thr.join();
    expected_thr.join();

    test_wls_pack(actual_pack, expected_pack);
}

INSTANTIATE_TEST_SUITE_P(
    WLSSuite, wls_fixture,

    // combination of inputs: (seed, n, p, intr, maxit)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(5, 10, 20, 30, 40),
        testing::Values(0, 1),
        testing::Values(1, 50, 100)
        ),
    
    // more informative test name
    [](const testing::TestParamInfo<wls_fixture::ParamType>& info) {
      std::string name = 
          std::string("seed_") + std::to_string(std::get<0>(info.param)) + '_'
          + "n_" + std::to_string(std::get<1>(info.param)) + '_'
          + "p_" + std::to_string(std::get<2>(info.param)) + '_'
          + "intr_" + std::to_string(std::get<3>(info.param)) + '_'
          + "maxit_" + std::to_string(std::get<4>(info.param));
      return name;
    }

);

} // namespace glmnetpp
