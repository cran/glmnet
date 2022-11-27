#include <benchmark/benchmark.h>
#include <testutil/data_util.hpp>
#include <testutil/mock_pb.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_driver/binomial.hpp>
#include <glmnetpp_bits/elnet_path.hpp>
#include <glmnetpp_bits/elnet_point.hpp>

namespace glmnetpp {

template <bool do_glmnetpp>
struct binomial_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::binomial>;

    elnet_driver_t elnet_driver;

    bool isd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::VectorXd ca;
    Eigen::MatrixXd X, X_cache, y, y_cache, g, g_cache, a0, cl, cl_cache;
    Eigen::VectorXd dev, alm;
    Eigen::VectorXi ia, nin; 
    double dev0 = 0;

    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(g_cache.cols(), nlam);
        ca.setZero(nx * g_cache.cols() * nlam);
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
        y = y_cache;
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
        int kopt = state.range(2);
        bool do_two_class = state.range(3);

        DataGen dgen(seed);
        X_cache = dgen.make_X(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p);
        y_cache.resize(n, 2);
        y_cache.col(0) = dgen.make_y(X_cache, beta);
        y_cache.col(0).array() = (y_cache.col(0).array() > y_cache.col(0).mean()).template cast<double>();
        y_cache.col(1).array() = 1. - y_cache.col(0).array();
        auto ne = dgen.make_ne(p);
        auto nx = p;

        if (do_two_class) {
            g_cache.setZero(n, 1);
        } else {
            g_cache.setZero(n, 2);
        }

        init(p, nx, nlam);

        state.counters["glmnetpp"] = do_glmnetpp;
        state.counters["kopt"] = kopt;
        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["sp"] = false;
        state.counters["two_class"] = do_two_class;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_driver.fit(
                        alpha, X, y, g, jd, vp, cl, ne, nx,
                        nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                        a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int no = X.rows();
                int ni = X.cols();
                int nc = g.cols();
                int iisd = isd;
                int iintr = intr;
                ::lognet_(
                        const_cast<double*>(&alpha), &no, &ni, &nc,
                        X.data(), y.data(), g.data(),
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
                        const_cast<int*>(&kopt), &lmu, a0.data(), ca.data(), ia.data(), 
                        nin.data(), &dev0, dev.data(), alm.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        binomial_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        binomial_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(binomial_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1, 2},
        {true}
        })  
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1},
        {false}
        })
    ;

BENCHMARK_REGISTER_F(binomial_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1, 2},
        {true}
        })  
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1},
        {false}
        })
    ;

// ==========================================================================

template <bool do_glmnetpp>
struct sp_binomial_fixture : benchmark::Fixture
{
    using elnet_driver_t = ElnetDriver<util::glm_type::binomial>;

    elnet_driver_t elnet_driver;

    bool isd = true, intr = true;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::VectorXd ca;
    Eigen::SparseMatrix<double> X, X_cache;
    Eigen::MatrixXd y, y_cache, g, g_cache, a0, cl, cl_cache;
    Eigen::VectorXd dev, alm;
    Eigen::VectorXi ia, nin, x_inner, x_outer;
    double dev0 = 0;

    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        a0.setZero(g_cache.cols(), nlam);
        ca.setZero(nx * g_cache.cols() * nlam);
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
        y = y_cache;
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
        int kopt = state.range(2);
        bool do_two_class = state.range(3);

        DataGen dgen(seed);
        X_cache = dgen.make_X_sparse(n, p);
        cl_cache = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        auto jd = dgen.make_jd(p);
        auto beta = dgen.make_beta(p);
        y_cache.resize(n, 2);
        y_cache.col(0) = dgen.make_y(X_cache, beta);
        y_cache.col(0).array() = (y_cache.col(0).array() > y_cache.col(0).mean()).template cast<double>();
        y_cache.col(1).array() = 1. - y_cache.col(0).array();
        auto ne = dgen.make_ne(p);
        auto nx = p;

        if (do_two_class) {
            g_cache.setZero(n, 1);
        } else {
            g_cache.setZero(n, 2);
        }

        init(p, nx, nlam);

        state.counters["glmnetpp"] = do_glmnetpp;
        state.counters["kopt"] = kopt;
        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["sp"] = true;
        state.counters["two_class"] = do_two_class;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_driver.fit(
                        alpha, X, y, g, jd, vp, cl, ne, nx,
                        nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                        a0, ca, ia, nin, dev0, dev, alm, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int no = X.rows();
                int ni = X.cols();
                int nc = g.cols();
                int iisd = isd;
                int iintr = intr;
                ::splognet_(
                    const_cast<double*>(&alpha), &no, &ni, &nc,
                    const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(), y.data(), g.data(),
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
                    const_cast<int*>(&kopt), &lmu, a0.data(), ca.data(), ia.data(), 
                    nin.data(), &dev0, dev.data(), alm.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        sp_binomial_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        sp_binomial_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(sp_binomial_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1, 2},
        {true}
        })  
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1},
        {false}
        })
    ;

BENCHMARK_REGISTER_F(sp_binomial_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1, 2},
        {true}
        })  
    ->ArgsProduct({
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2),
        {0, 1},
        {false}
        })
    ;

} // namespace glmnetpp
