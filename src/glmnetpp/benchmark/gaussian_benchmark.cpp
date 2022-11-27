#include <benchmark/benchmark.h>
#include <testutil/data_util.hpp>
#include <testutil/mock_pb.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_driver/gaussian.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_cov.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/sp_gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_multi.hpp>

namespace glmnetpp {

/*
 * Dense naive/cov
 */
template <bool do_glmnetpp>
struct gaussian_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::gaussian>;

    elnet_driver_t elnet_driver;

    bool isd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::VectorXd y, y_cache, w, w_cache, a0;
    Eigen::MatrixXd X, X_cache, ca, cl, cl_cache;
    Eigen::VectorXd rsq, alm;
    Eigen::VectorXi ia, nin; 


    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(nlam);
        ca.setZero(nx, nlam);
        ia.setZero(p);
        nin.setZero(nlam);
        rsq.setZero(nlam);
        alm.setZero(nlam);
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        X = X_cache;
        y = y_cache;
        w = w_cache;
        cl = cl_cache;
        a0.setZero();
        ca.setZero();
        ia.setZero();
        nin.setZero();
        rsq.setZero();
        alm.setZero();
    }

    void run(benchmark::State& state)
    {
        int seed = 123124;
        int n = state.range(0);
        int p = state.range(1);
        bool ka = state.range(2);

        DataGen dgen(seed);
        X_cache = dgen.make_X(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        w_cache = dgen.make_w(n);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p);
        y_cache = dgen.make_y(X_cache, beta);
        auto ne = dgen.make_ne(p);
        auto nx = dgen.make_nx(p);

        init(p, nx, nlam);

        state.counters["glmnetpp"] = do_glmnetpp;
        state.counters["ka"] = ka;
        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["sp"] = false;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_driver.fit(ka, alpha, X, y, w, jd, vp, cl, ne, nx, 
                               nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                               a0, ca, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int ika = ka + 1;
                int ni = X.cols();
                int no = X.rows();
                int iisd = isd;
                int iintr = intr;
                ::elnet_(
                    &ika, 
                    &alpha, &no, &ni,
                    X.data(), y.data(), w.data(),
                    jd.data(), vp.data(), cl.data(), &ne, 
                    &nx, &nlam, &flmin, ulam.data(), &thr, &iisd, &iintr, &maxit, 
                    &lmu, a0.data(), ca.data(), ia.data(), 
                    nin.data(), rsq.data(), alm.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        gaussian_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        gaussian_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(gaussian_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {true, false}
        })
    ;

BENCHMARK_REGISTER_F(gaussian_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {true, false}
        })
    ;

// ==========================================================================

/*
 * Sparse naive/cov
 */
template <bool do_glmnetpp>
struct sp_gaussian_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::gaussian>;

    elnet_driver_t elnet_driver;

    bool isd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::VectorXd y, y_cache, w, w_cache, a0;
    Eigen::SparseMatrix<double> X, X_cache;
    Eigen::MatrixXd ca, cl, cl_cache;
    Eigen::VectorXd rsq, alm;
    Eigen::VectorXi ia, nin, x_inner, x_outer;

    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(nlam);
        ca.setZero(nx, nlam);
        ia.setZero(p);
        nin.setZero(nlam);
        rsq.setZero(nlam);
        alm.setZero(nlam);
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        X = X_cache;
        X.makeCompressed();
        x_inner = make_sp_inner_idx_1idx(X);
        x_outer = make_sp_outer_idx_1idx(X);
        y = y_cache;
        w = w_cache;
        cl = cl_cache;
        a0.setZero();
        ca.setZero();
        ia.setZero();
        nin.setZero();
        rsq.setZero();
        alm.setZero();
    }

    void run(benchmark::State& state)
    {
        int seed = 123124;
        int n = state.range(0);
        int p = state.range(1);
        bool ka = state.range(2);

        DataGen dgen(seed);
        X_cache = dgen.make_X_sparse(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        w_cache = dgen.make_w(n);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p);
        y_cache = dgen.make_y(X_cache, beta);
        auto ne = dgen.make_ne(p);
        auto nx = dgen.make_nx(p);

        init(p, nx, nlam);

        state.counters["glmnetpp"] = do_glmnetpp;
        state.counters["ka"] = ka;
        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["sp"] = true;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_driver.fit(ka, alpha, X, y, w, jd, vp, cl, ne, nx, 
                               nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                               a0, ca, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int ika = ka + 1;
                int ni = X.cols();
                int no = X.rows();
                int iisd = isd;
                int iintr = intr;
                ::spelnet_(
                    &ika, 
                    &alpha, &no, &ni,
                    X.valuePtr(), x_outer.data(), x_inner.data(),
                    y.data(), w.data(),
                    jd.data(), vp.data(), cl.data(), &ne, 
                    &nx, &nlam, &flmin, ulam.data(), &thr, &iisd, &iintr, &maxit, 
                    &lmu, a0.data(), ca.data(), ia.data(), 
                    nin.data(), rsq.data(), alm.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        sp_gaussian_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        sp_gaussian_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(sp_gaussian_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {true, false}
        })
    ;

BENCHMARK_REGISTER_F(sp_gaussian_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {true, false}
        })
    ;

// ==========================================================================

/*
 * Dense multi
 */
template <bool do_glmnetpp>
struct gaussian_multi_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::gaussian>;

    elnet_driver_t elnet_driver;

    bool isd = true, jsd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::MatrixXd y, y_cache, a0;
    Eigen::VectorXd w, w_cache, ca;
    Eigen::MatrixXd X, X_cache, cl, cl_cache;
    Eigen::VectorXd rsq, alm;
    Eigen::VectorXi ia, nin; 

    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(y_cache.cols(), nlam);
        ca.setZero(nx * y_cache.cols() * nlam);
        ia.setZero(p);
        nin.setZero(nlam);
        rsq.setZero(nlam);
        alm.setZero(nlam);
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        X = X_cache;
        y = y_cache;
        w = w_cache;
        cl = cl_cache;
        a0.setZero();
        ca.setZero();
        ia.setZero();
        nin.setZero();
        rsq.setZero();
        alm.setZero();
    }

    void run(benchmark::State& state)
    {
        int seed = 123124;
        int n = state.range(0);
        int p = state.range(1);
        int nr = 3;

        DataGen dgen(seed);
        X_cache = dgen.make_X(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        w_cache = dgen.make_w(n);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p);
        y_cache.resize(n, nr);
        for (int i = 0; i < nr; ++i) {
            y_cache.col(i) = dgen.make_y(X_cache, beta);
        }
        auto ne = dgen.make_ne(p);
        auto nx = dgen.make_nx(p);

        init(p, nx, nlam);

        state.counters["glmnetpp"] = do_glmnetpp;
        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["sp"] = false;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_driver.fit(
                        alpha, X, y, w, jd, vp, cl, ne, nx, 
                        nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                        a0, ca, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int ni = X.cols();
                int nr = y.cols();
                int no = X.rows();
                int iisd = isd;
                int ijsd = jsd;
                int iintr = intr;
                multelnet_(
                        const_cast<double*>(&alpha), &no, &ni, &nr,
                        X.data(), 
                        y.data(),
                        w.data(),
                        const_cast<int*>(jd.data()), 
                        const_cast<double*>(vp.data()), 
                        const_cast<double*>(cl.data()),
                        const_cast<int*>(&ne), 
                        const_cast<int*>(&nx), 
                        const_cast<int*>(&nlam),
                        const_cast<double*>(&flmin), 
                        const_cast<double*>(ulam.data()), 
                        const_cast<double*>(&thr), 
                        &iisd, &ijsd, &iintr,
                        const_cast<int*>(&maxit), 
                        &lmu,
                        a0.data(), ca.data(), ia.data(), nin.data(), rsq.data(), alm.data(),
                        &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        gaussian_multi_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        gaussian_multi_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(gaussian_multi_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

BENCHMARK_REGISTER_F(gaussian_multi_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

// ==========================================================================

/*
 * Sparse multi
 */
template <bool do_glmnetpp>
struct sp_gaussian_multi_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::gaussian>;

    elnet_driver_t elnet_driver;

    bool isd = true, jsd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::MatrixXd y, y_cache, a0;
    Eigen::VectorXd w, w_cache, ca;
    Eigen::SparseMatrix<double> X, X_cache;
    Eigen::MatrixXd cl, cl_cache;
    Eigen::VectorXd rsq, alm;
    Eigen::VectorXi ia, nin, x_inner, x_outer;


    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(y_cache.cols(), nlam);
        ca.setZero(nx * y_cache.cols() * nlam);
        ia.setZero(p);
        nin.setZero(nlam);
        rsq.setZero(nlam);
        alm.setZero(nlam);
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        X = X_cache;
        X.makeCompressed();
        x_inner = make_sp_inner_idx_1idx(X);
        x_outer = make_sp_outer_idx_1idx(X);
        y = y_cache;
        w = w_cache;
        cl = cl_cache;
        a0.setZero();
        ca.setZero();
        ia.setZero();
        nin.setZero();
        rsq.setZero();
        alm.setZero();
    }

    void run(benchmark::State& state)
    {
        int seed = 123124;
        int n = state.range(0);
        int p = state.range(1);
        int nr = 3;

        DataGen dgen(seed);
        X_cache = dgen.make_X_sparse(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        w_cache = dgen.make_w(n);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p);
        y_cache.resize(n, nr);
        for (int i = 0; i < nr; ++i) {
            y_cache.col(i) = dgen.make_y(X_cache, beta);
        }
        auto ne = dgen.make_ne(p);
        auto nx = dgen.make_nx(p);

        init(p, nx, nlam);

        state.counters["glmnetpp"] = do_glmnetpp;
        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["sp"] = true;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_driver.fit(
                        alpha, X, y, w, jd, vp, cl, ne, nx, 
                        nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                        a0, ca, ia, nin, rsq, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int ni = X.cols();
                int nr = y.cols();
                int no = X.rows();
                int iisd = isd;
                int ijsd = jsd;
                int iintr = intr;
                multspelnet_(
                        const_cast<double*>(&alpha), &no, &ni, &nr,
                        const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(), 
                        y.data(),
                        w.data(),
                        const_cast<int*>(jd.data()), 
                        const_cast<double*>(vp.data()), 
                        const_cast<double*>(cl.data()),
                        const_cast<int*>(&ne), 
                        const_cast<int*>(&nx), 
                        const_cast<int*>(&nlam),
                        const_cast<double*>(&flmin), 
                        const_cast<double*>(ulam.data()), 
                        const_cast<double*>(&thr), 
                        &iisd, &ijsd, &iintr,
                        const_cast<int*>(&maxit), 
                        &lmu,
                        a0.data(), ca.data(), ia.data(), nin.data(), rsq.data(), alm.data(),
                        &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        sp_gaussian_multi_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        sp_gaussian_multi_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(sp_gaussian_multi_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

BENCHMARK_REGISTER_F(sp_gaussian_multi_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

} // namespace glmnetpp
