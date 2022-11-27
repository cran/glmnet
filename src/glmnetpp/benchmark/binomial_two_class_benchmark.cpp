#include <benchmark/benchmark.h>
#include <testutil/data_util.hpp>
#include <testutil/mock_pb.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_point/binomial_two_class.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp>
#include <glmnetpp_bits/elnet_point/internal/binomial_two_class.hpp>

namespace glmnetpp {

template <bool do_glmnetpp>
struct binomial_two_class_fixture : benchmark::Fixture
{
    using internal_t = ElnetPointInternal<
                util::glm_type::binomial,
                util::mode_type<util::glm_type::binomial>::two_class,
                double,
                int,
                int>;
    using elnet_point_t = ElnetPoint<
            util::glm_type::binomial,
            util::mode_type<util::glm_type::binomial>::two_class,
            internal_t>;
    using elnet_path_t = ElnetPath<
        util::glm_type::binomial,
        util::mode_type<util::glm_type::binomial>::two_class,
        elnet_point_t>;

    elnet_path_t elnet_path;

    Eigen::VectorXd xm, xs;
    double flmin = 0.;
    double alpha = 1.0, thr = 1e-14;
    int ne, maxit = 100000, nlam = 100, kopt = 0; // kopt doesn't matter for two-class 
    bool isd = 1, intr = 1;
    Eigen::VectorXd ulam;

    int lmu, nlp, jerr;
    double dev0;
    Eigen::MatrixXd ao;
    Eigen::VectorXd a0, almo, dev, g;
    Eigen::VectorXi ia, kin; 

    void init(int p, int n, int nx, int nlam)
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        dev0 = 0;
        ao.setZero(nx, nlam);
        a0.setZero(nlam);
        ia.setZero(p);
        kin.setZero(nlam);
        dev.setZero(nlam);
        almo.setZero(nlam);
        g.setZero(n);
    }

    void reset()
    {
        lmu = 0;
        nlp = 0;
        jerr = 0;
        dev0 = 0;
        ao.setZero();
        a0.setZero();
        ia.setZero();
        kin.setZero();
        dev.setZero();
        almo.setZero();
        g.setZero();
    }

    void run(benchmark::State& state)
    {
        int seed = 1424821;
        int n = state.range(0);
        int p = state.range(1);

        DataGen dgen(seed);
        auto X = dgen.make_X(n, p);
        auto cl = dgen.make_cl(p);
        auto vp = dgen.make_vp(p);
        Eigen::VectorXd w(X.rows()); w.array() = 1./w.size();
        auto ju = dgen.make_ju(p);
        auto beta = dgen.make_beta(p);
        Eigen::VectorXd y = dgen.make_y(X, beta);
        y.array() =  (y.array() > 0.).template cast<double>();
        Eigen::VectorXd g(n); g.setZero();
        auto nx = dgen.make_nx(p);
        Eigen::VectorXi ia(nx); ia.setZero();

        vp /= vp.sum() / p;
        Chkvars::eval(X, ju);

        xm.setZero(p);
        xs.setZero(p);
        init(p, n, nx, nlam);

        LStandardize1::eval(X, w, ju, isd, intr, xm, xs);

        for (int j = 0; j < p; ++j) {
            cl.col(j) *= xs(j);
        }
        ne = p;

        state.counters["n"] = n;
        state.counters["p"] = p;
        state.counters["glmnetpp"] = do_glmnetpp;

        //if constexpr (do_glmnetpp) {
        //    reset();
        //    Eigen::MatrixXd ao_copy_1 = ao, ao_copy_2 = ao;
        //    int lmu_copy_1 = 0, lmu_copy_2=0, nlp_copy_1=0, nlp_copy_2=0, jerr_copy_1=0, jerr_copy_2=0;
        //    double dev0_copy_1=0, dev0_copy_2=0;
        //    Eigen::VectorXd g_copy_1=g, g_copy_2=g, a0_copy_1=a0, a0_copy_2=a0, almo_copy_1=almo, almo_copy_2=almo, dev_copy_1=dev, dev_copy_2=dev;
        //    Eigen::VectorXi ia_copy_1=ia, ia_copy_2=ia, kin_copy_1=kin, kin_copy_2=kin; 

        //    elnet_path.fit(alpha, ju, vp, cl, ne, nx, X, y, g_copy_1, w, 
        //                   nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu_copy_1,
        //                   a0_copy_1, ao_copy_1, ia_copy_1, kin_copy_1, dev0_copy_1, dev_copy_1, almo_copy_1, 
        //                   nlp_copy_1, jerr_copy_1, mock_setpb, InternalParams());
        //    int ni = X.cols();
        //    int no = X.rows();
        //    int iisd = isd;
        //    int iintr = intr;
        //    lognet2n_(&alpha, &no, &ni, X.data(), y.data(), g_copy_2.data(), w.data(),
        //            ju.data(), vp.data(), cl.data(),
        //            &ne, &nx, &nlam,
        //            &flmin, ulam.data(), &thr, &iisd, &iintr, &maxit, &kopt, 
        //            &lmu_copy_2, a0_copy_2.data(), ao_copy_2.data(), ia_copy_2.data(), kin_copy_2.data(), &dev0_copy_2, dev_copy_2.data(),
        //            almo_copy_2.data(), &nlp_copy_2, &jerr_copy_2);
        //    std::cout << std::abs(lmu_copy_1 - lmu_copy_2) << ',';
        //    std::cout << ((a0_copy_1.array() - a0_copy_2.array()).abs()).maxCoeff() << ',';
        //    std::cout << ((ao_copy_1.array() - ao_copy_2.array()).abs()).maxCoeff() << ',';
        //    std::cout << ((ia_copy_1.array() - ia_copy_2.array()).abs()).maxCoeff() << ',';
        //    std::cout << ((kin_copy_1.array() - kin_copy_2.array()).abs()).maxCoeff() << ',';
        //    std::cout << std::abs(dev0_copy_1 - dev0_copy_2) << ',';
        //    std::cout << ((dev_copy_1.array() - dev_copy_2.array()).abs()).maxCoeff() << ',';
        //    std::cout << ((almo_copy_1.array() - almo_copy_2.array()).abs()).maxCoeff() << ',';
        //    std::cout << std::abs(nlp_copy_1 - nlp_copy_2) << ',';
        //    std::cout << nlp_copy_1 << ',';
        //    std::cout << std::abs(jerr_copy_1 - jerr_copy_2) << std::endl;
        //    //std::cout << "glmnetpp: " << n << ',' << p << '\n' << ao_copy_1 << std::endl;
        //    //std::cout << "glmnet: " << n << ',' << p << '\n' << ao_copy_2 << std::endl;
        //}

        for (auto _ : state) {
            state.PauseTiming();
            reset();
            state.ResumeTiming();
            if constexpr (do_glmnetpp) {
                elnet_path.fit(alpha, ju, vp, cl, ne, nx, X, y, g, w, 
                               nlam, flmin, ulam, thr, isd, intr, maxit, kopt, lmu,
                               a0, ao, ia, kin, dev0, dev, almo, 
                               nlp, jerr, mock_setpb, InternalParams());
            } else {
                int ni = X.cols();
                int no = X.rows();
                int iisd = isd;
                int iintr = intr;
                lognet2n_(&alpha, &no, &ni, X.data(), y.data(), g.data(), w.data(),
                        ju.data(), vp.data(), cl.data(),
                        &ne, &nx, &nlam,
                        &flmin, ulam.data(), &thr, &iisd, &iintr, &maxit, &kopt, 
                        &lmu, a0.data(), ao.data(), ia.data(), kin.data(), &dev0, dev.data(),
                        almo.data(), &nlp, &jerr);
            }
        }
    }
};


BENCHMARK_TEMPLATE_DEFINE_F(
        binomial_two_class_fixture,
        glmnetpp,
        true)(benchmark::State& state)
{ run(state); }

BENCHMARK_TEMPLATE_DEFINE_F(
        binomial_two_class_fixture,
        legacy,
        false)(benchmark::State& state)
{ run(state); }

BENCHMARK_REGISTER_F(binomial_two_class_fixture,
                     glmnetpp)
    ->ArgsProduct({
        {1000},
        //{1024}
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

BENCHMARK_REGISTER_F(binomial_two_class_fixture,
                     legacy)
    ->ArgsProduct({
        {1000},
        //{1024}
        benchmark::CreateRange(1<<5, 1<<11, 2)
        })
    ;

} // namespace glmnetpp
