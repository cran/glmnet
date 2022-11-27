#include <benchmark/benchmark.h>
#include <testutil/data_util.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/internal.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/gaussian_naive.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_naive.hpp>

namespace glmnetpp {

template <bool do_glmnetpp>
struct gaussian_naive_fixture : benchmark::Fixture
{
    using internal_t = ElnetPointInternal<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::naive,
                double,
                int,
                int>;
    using elnet_point_t = ElnetPoint<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::naive,
            internal_t>;
    using elnet_path_t = ElnetPath<
        util::glm_type::gaussian,
        util::mode_type<util::glm_type::gaussian>::naive,
        elnet_point_t>;

    elnet_path_t elnet_path;

    Eigen::VectorXd xm, xs, xv;
    double ym = 0, ys = 0, flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int ne, maxit = 100000, nlam = 100;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    Eigen::MatrixXd ao;
    Eigen::VectorXd rsqo, almo;
    Eigen::VectorXi ia, kin; 

    void init(int p, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        ao.setZero(nx, nlam);
        ia.setZero(p);
        kin.setZero(nlam);
        rsqo.setZero(nlam);
        almo.setZero(nlam);
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        ao.setZero();
        ia.setZero();
        kin.setZero();
        rsqo.setZero();
        almo.setZero();
    }

    void run(benchmark::State& state)
    {
        int seed = 123124;
        int n = state.range(0);
        int p = state.range(1);

        DataGen dgen(seed);
        auto X = dgen.make_X(n, p);
        auto cl = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        auto w = dgen.make_w(n);
        auto ju = dgen.make_ju(p);
        auto beta = dgen.make_beta(p);
        auto y = dgen.make_y(X, beta);
        Eigen::VectorXd g(p); g.setZero();
        Eigen::VectorXd xv(p); xv.setOnes();
        auto nx = dgen.make_nx(p);
        Eigen::VectorXi ia(nx); ia.setZero();

        vp /= vp.sum() / p;
        Chkvars::eval(X, ju);

        xm.setZero(p);
        xs.setZero(p);
        xv.setZero(p);
        init(p, nx, nlam);

        Standardize1::eval(X, y, w, false, true, 
                    ju, xm, xs, ym, ys, xv);
        cl /= ys; 
        for (int j = 0; j < p; ++j) {
            cl.col(j) *= xs(j);
        }
        ne = p;

        state.counters["n"] = n;
        state.counters["p"] = p;

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_path.fit(alpha, ju, vp, cl, y, ne, nx, X,
                               nlam, flmin, ulam, thr, maxit, xv, lmu,
                               ao, ia, kin, rsqo, almo, nlp, jerr, mock_setpb, InternalParams());
            } else {
                int ni = X.cols();
                int no = X.rows();
                elnet2_(&alpha, &ni, ju.data(), vp.data(), cl.data(),
                        y.data(), &no, &ne, &nx, X.data(), &nlam,
                        &flmin, ulam.data(), &thr, &maxit, xv.data(), &lmu,
                        ao.data(), ia.data(), kin.data(), rsqo.data(), almo.data(),
                        &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        gaussian_naive_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        gaussian_naive_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(gaussian_naive_fixture,
                     glmnetpp)
    ->ArgsProduct({
        //{100, 500, 1000, 2000},
        //benchmark::CreateRange(2, 1<<11, 2)
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

BENCHMARK_REGISTER_F(gaussian_naive_fixture,
                     legacy)
    ->ArgsProduct({
        //{100, 500, 1000, 2000},
        //benchmark::CreateRange(2, 1<<11, 2)
        {1000},
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

} // namespace glmnetpp
