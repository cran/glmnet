#include <elnet_path/gaussian_base.hpp>
#include <testutil/mock_pb.hpp>
#include <testutil/translation/multspelnet2.hpp>
#include <testutil/internal.hpp>
#include <glmnetpp_bits/elnet_point/internal/sp_gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_point/internal/gaussian_base.hpp> // naive
#include <glmnetpp_bits/elnet_point/sp_gaussian_multi.hpp>
#include <glmnetpp_bits/elnet_path/base.hpp>
#include <glmnetpp_bits/elnet_path/sp_gaussian_multi.hpp>

namespace glmnetpp {

struct SpGaussianMultiPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;

    // will be modified
    Eigen::MatrixXd y;
    Eigen::VectorXd ao;
    Eigen::VectorXi ia;
    Eigen::VectorXi kin;
    Eigen::VectorXd rsqo; 
    Eigen::VectorXd almo;
    double ys0;
    int nlp = 0, jerr = 0, lmu = 0;
    
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& w;
    const Eigen::VectorXd& xm;
    const Eigen::VectorXd& xs;
    const Eigen::VectorXd& xv;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXd& cl;
    const Eigen::VectorXi& ju;

    SpGaussianMultiPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            double _alpha,
            double _flmin,
            double _ys0,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _xm,
            const Eigen::VectorXd& _xs,
            const Eigen::VectorXd& _xv,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::VectorXd& _cl,
            const Eigen::VectorXi& _ju)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha), flmin(_flmin)
        , y(_y)
        , ao(_nx * _y.cols() * _nlam)
        , ia(nx)
        , kin(nlam)
        , rsqo(nlam)
        , almo(nlam)
        , ys0(_ys0)
        , X(_X), w(_w), xm(_xm), xs(_xs), xv(_xv), ulam(_ulam), vp(_vp), cl(_cl), ju(_ju)
    {
        ao.setZero();
        ia.setZero();
        kin.setZero();
        rsqo.setZero();
        almo.setZero();
    }

    void fit() 
    {
        using internal_t = SpElnetPointInternal<
                    util::glm_type::gaussian,
                    util::mode_type<util::glm_type::gaussian>::multi,
                    double,
                    int, 
                    int>;
        using elnet_point_t = SpElnetPoint<
                util::glm_type::gaussian,
                util::mode_type<util::glm_type::gaussian>::multi,
                internal_t>;
        using elnet_path_t = SpElnetPath<
            util::glm_type::gaussian,
            util::mode_type<util::glm_type::gaussian>::multi,
            elnet_point_t>;

        elnet_path_t elnet_path;
        elnet_path.fit(
                alpha, ju, vp, cl, y, w, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, ys0, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr, mock_setpb, InternalParams());
    }

    void fit_old() 
    {
        transl::multspelnet2<double>(
                alpha, ju, vp, cl, y, w, ne, nx, X,
                nlam, flmin, ulam, thr, maxit, xm, xs, xv, ys0, lmu,
                ao, ia, kin, rsqo, almo, nlp, jerr);
    }
};

struct sp_gaussian_multi_fixture
    : gaussian_fixture
{
    void SetUp() override
    {
        size_t seed, n, p, nr = 3;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p);
        y.resize(n, nr);
        for (size_t i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl.resize(2 * nr * p);
        auto cl_sub = dgen.make_cl(p);
        for (size_t i = 0; i < p; ++i) {
            Eigen::Map<Eigen::MatrixXd> cl_slice(
                    cl.data() + i * 2 * nr, 2, nr);
            cl_slice.colwise() = cl_sub.col(i);
        }
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        xv.setOnes(p);

        xm.setZero(p);
        xs.setZero(p);
        Eigen::VectorXd ym(nr), ys(nr);
        MultSpStandardize1::eval(X, y, w, 1, 1, 1, ju, xm, xs, ym, ys, xv, ys0);
    }

protected:
    using base_t = gaussian_fixture;
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd y;
    Eigen::VectorXd cl, w, xm, xs;
    double ys0;

    void check_pack(const SpGaussianMultiPack& actual,
                    const SpGaussianMultiPack& expected)
    {
        expect_near_vec(actual.ao, expected.ao, 1e-15);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.kin, expected.kin);
        expect_double_eq_vec(actual.rsqo, expected.rsqo);
        expect_double_eq_vec(actual.almo, expected.almo);
        EXPECT_DOUBLE_EQ(actual.ys0, expected.ys0);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
        expect_double_eq_mat(actual.y, expected.y);
    }
};

TEST_P(sp_gaussian_multi_fixture, sp_gaussian_multi_test)
{
    SpGaussianMultiPack actual(
            maxit, nx, ne, nlam, alpha, flmin, ys0,
            X, y, w, xm, xs, xv, ulam, vp, cl, ju);
    SpGaussianMultiPack expected(actual);
    run(actual, expected, 4, 5); 
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpGaussianMultiSuite, sp_gaussian_multi_fixture,
    testing::Combine(
        testing::Values(241, 412, 23968, 31, 87201, 746, 104),
        testing::Values(10, 30, 50),
        testing::Values(5, 20, 40, 60),
        testing::Values(1, 50, 100),
        testing::Values(1, 4, 10),
        testing::Values(0.0, 0.5, 1.0),
        testing::Values(0.5, 1.0, 1.5)
        )
);

} // namespace glmnetpp
