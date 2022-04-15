#include <legacy/legacy.h>
#include <testutil/translation/multspelnet.hpp>
#include <translation/elnet_base_fixture.hpp>

namespace glmnetpp {

struct MultSpElnetPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;
    const bool isd, jsd, intr;

    // will be modified
    Eigen::MatrixXd y;
    Eigen::VectorXd w;
    Eigen::MatrixXd a0;
    Eigen::VectorXd ao;
    Eigen::VectorXi ia;
    Eigen::VectorXi kin;
    Eigen::VectorXd rsqo; 
    Eigen::VectorXd almo;
    int nlp = 0, jerr = 0, lmu = 0;
    
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::VectorXi& jd;
    const Eigen::MatrixXd& cl;

    MultSpElnetPack(
            int _maxit,
            int _nx,
            int _ne,
            int _nlam,
            bool _isd,
            bool _jsd,
            bool _intr,
            double _alpha,
            double _flmin,
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXd& _ulam,
            const Eigen::VectorXd& _vp,
            const Eigen::MatrixXd& _cl,
            const Eigen::VectorXi& _jd)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam), alpha(_alpha), flmin(_flmin), isd(_isd), jsd(_jsd), intr(_intr)
        , y(_y)
        , w(_w)
        , a0(_y.cols(), _nlam)
        , ao(_nx * _y.cols() * _nlam)
        , ia(nx)
        , kin(nlam)
        , rsqo(nlam)
        , almo(nlam)
        , X(_X)
        , ulam(_ulam), vp(_vp), jd(_jd)
        , cl(_cl)
    {
        ao.setZero();
        ia.setZero();
        kin.setZero();
        rsqo.setZero();
        almo.setZero();
    }

    void fit()
    {
        transl::multspelnet<float>(
                alpha, X, y, w, jd, vp, cl, ne, nx, 
                nlam, flmin, ulam, thr, isd, jsd, intr, maxit, lmu,
                a0, ao, ia, kin, rsqo, almo, nlp, jerr);
    }

    void fit_old() 
    {
        int ni = X.cols();
        int nr = y.cols();
        int no = X.rows();
        int iisd = isd;
        int ijsd = jsd;
        int iintr = intr;
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
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
                a0.data(), ao.data(), ia.data(), kin.data(), rsqo.data(), almo.data(),
                &nlp, &jerr);
    }
};

struct multspelnet_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double, bool, bool, bool> >
{
    void SetUp() override
    {
        size_t seed, n, p, nr = 3;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin, isd, jsd, intr) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p);
        y.resize(n, nr);
        for (size_t i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        w = dgen.make_w(n);
        jd = dgen.make_jd(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
    }

protected:
    using base_t = elnet_base_fixture;
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd y, cl;
    Eigen::VectorXd w, vp, ulam;
    Eigen::VectorXi jd;
    size_t maxit, nlam, nx, ne;
    double alpha, flmin;
    bool isd, jsd, intr;

    void check_pack(const MultSpElnetPack& actual,
                    const MultSpElnetPack& expected)
    {
        expect_near_mat(actual.y, expected.y, 1e-7);
        expect_float_eq_vec(actual.w, expected.w);
        expect_near_mat(actual.a0, expected.a0, 1e-7);
        expect_near_vec(actual.ao, expected.ao, 1e-7);
        expect_eq_vec(actual.ia, expected.ia);
        expect_eq_vec(actual.kin, expected.kin);
        expect_near_vec(actual.rsqo, expected.rsqo, 1e-7);
        expect_float_eq_vec(actual.almo, expected.almo);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);
        EXPECT_EQ(actual.lmu, expected.lmu);
    }
};

TEST_P(multspelnet_fixture, multspelnet_test)
{
    MultSpElnetPack actual(
            maxit, nx, ne, nlam, isd, jsd, intr, alpha, flmin,
            X, y, w, ulam, vp, cl, jd);
    MultSpElnetPack expected(actual);
    run(actual, expected, 4, 5);
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultSpElnetSuite, multspelnet_fixture,
    testing::Combine(
        testing::Values(241, 412, 3, 31), // seed
        testing::Values(10, 30, 50),          // n
        testing::Values(5, 20, 40),           // p
        testing::Values(1, 50),               // maxit
        testing::Values(1, 4),                // nlam
        testing::Values(0.0, 0.5, 1.0),       // alpha
        testing::Values(0.5, 1.0, 1.5),       // flmin
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

} // namespace glmnetpp
