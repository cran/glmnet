#include <benchmark/benchmark.h>
#include <testutil/data_util.hpp>
#include <testutil/mock_pb.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_driver/poisson.hpp>
#include <glmnetpp_bits/elnet_path.hpp>
#include <glmnetpp_bits/elnet_point.hpp>

namespace glmnetpp {

template <bool do_glmnetpp>
struct poisson_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::poisson>;

    elnet_driver_t elnet_driver;

    bool isd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::VectorXd a0, y, g, g_cache;
    Eigen::MatrixXd X, X_cache, ca, cl, cl_cache;
    Eigen::VectorXd dev, alm;
    Eigen::VectorXi ia, nin; 
    double dev0 = 0;
    
    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(nlam);
        ca.setZero(nx, nlam);
        ia.setZero(p);
        nin.setZero(nlam);
        dev.setZero(nlam);
        alm.setZero(nlam);
        dev0 = 0;
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        X = X_cache;
        g = g_cache;
        cl = cl_cache;
        a0.setZero();
        ca.setZero();
        ia.setZero();
        nin.setZero();
        dev.setZero();
        alm.setZero();
        dev0 = 0;
    }

    void run(benchmark::State& state)
    {
        int seed = 1424821;
        int n = state.range(0);
        int p = state.range(1);

        DataGen dgen(seed);
        X_cache = dgen.make_X(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p, 0.01);
        Eigen::VectorXd f = (X_cache * beta).array().exp().matrix();
        std::mt19937 gen(seed);
        y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
        auto w = dgen.make_w(n);
        auto ne = dgen.make_ne(p);
        auto nx = p;
        g_cache.setZero(n);

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
                        alpha, X, y, g, w, jd, vp, cl, ne, nx,
                        nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                        a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int no = X.rows();
                int ni = X.cols();
                int iisd = isd;
                int iintr = intr;
                ::fishnet_(
                        const_cast<double*>(&alpha), &no, &ni,
                        X.data(), 
                        const_cast<double*>(y.data()), g.data(), 
                        const_cast<double*>(w.data()),
                        const_cast<int*>(jd.data()),
                        const_cast<double*>(vp.data()),
                        cl.data(),
                        const_cast<int*>(&ne), 
                        const_cast<int*>(&nx), 
                        const_cast<int*>(&nlam), 
                        const_cast<double*>(&flmin), 
                        const_cast<double*>(ulam.data()), 
                        const_cast<double*>(&thr), &iisd, &iintr, 
                        const_cast<int*>(&maxit), 
                        &lmu, a0.data(), ca.data(), ia.data(), 
                        nin.data(), &dev0, dev.data(), alm.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        poisson_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        poisson_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(poisson_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        })  
    ;

BENCHMARK_REGISTER_F(poisson_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })  
    ;

// ==========================================================================

template <bool do_glmnetpp>
struct sp_poisson_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::poisson>;

    elnet_driver_t elnet_driver;

    bool isd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::SparseMatrix<double> X, X_cache;
    Eigen::VectorXd a0, y, g, g_cache;
    Eigen::MatrixXd ca, cl, cl_cache;
    Eigen::VectorXd dev, alm;
    Eigen::VectorXi ia, nin, x_inner, x_outer; 
    double dev0 = 0;

    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(nlam);
        ca.setZero(nx, nlam);
        ia.setZero(p);
        nin.setZero(nlam);
        dev.setZero(nlam);
        alm.setZero(nlam);
        dev0 = 0;
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
        g = g_cache;
        cl = cl_cache;
        a0.setZero();
        ca.setZero();
        ia.setZero();
        nin.setZero();
        dev.setZero();
        alm.setZero();
        dev0 = 0;
    }

    void run(benchmark::State& state)
    {
        int seed = 1424821;
        int n = state.range(0);
        int p = state.range(1);

        DataGen dgen(seed);
        X_cache = dgen.make_X_sparse(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p, 0.01);
        Eigen::VectorXd f = (X_cache * beta).array().exp().matrix();
        std::mt19937 gen(seed);
        y = f.unaryExpr([&](auto m) -> double { return std::poisson_distribution<int>(m)(gen); });
        auto w = dgen.make_w(n);
        auto ne = dgen.make_ne(p);
        auto nx = p;
        g_cache.setZero(n);

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
                        alpha, X, y, g, w, jd, vp, cl, ne, nx,
                        nlam, flmin, ulam, thr, isd, intr, maxit, lmu,
                        a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int no = X.rows();
                int ni = X.cols();
                int iisd = isd;
                int iintr = intr;
                ::spfishnet_(
                        const_cast<double*>(&alpha), &no, &ni,
                        const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(), 
                        const_cast<double*>(y.data()), g.data(), 
                        const_cast<double*>(w.data()),
                        const_cast<int*>(jd.data()),
                        const_cast<double*>(vp.data()),
                        cl.data(),
                        const_cast<int*>(&ne), 
                        const_cast<int*>(&nx), 
                        const_cast<int*>(&nlam), 
                        const_cast<double*>(&flmin), 
                        const_cast<double*>(ulam.data()), 
                        const_cast<double*>(&thr), &iisd, &iintr, 
                        const_cast<int*>(&maxit), 
                        &lmu, a0.data(), ca.data(), ia.data(), 
                        nin.data(), &dev0, dev.data(), alm.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        sp_poisson_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        sp_poisson_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(sp_poisson_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })  
    ;

BENCHMARK_REGISTER_F(sp_poisson_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })  
    ;

} // namespace glmnetpp
