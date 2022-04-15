#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>
#include <tuple>

namespace glmnetpp {

// ======================================================
// Elnet Standardize
// ======================================================

struct Standardize1Pack
{
    Eigen::MatrixXd X;
    Eigen::VectorXd y, w, xm, xs, xv;
    Eigen::VectorXi ju;
    bool isd, intr;
    double ym = 0, ys = 0;
    int jerr = 0;

    Standardize1Pack(const Eigen::MatrixXd& _X,
                    const Eigen::VectorXd& _y,
                    const Eigen::VectorXd& _w,
                    const Eigen::VectorXi& _ju,
                    bool _isd, bool _intr)
        : X(_X), y(_y), w(_w), xm(_ju.size())
        , xs(_ju.size()), xv(_ju.size()), ju(_ju), isd(_isd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
        xv.setZero();
    }

    void standardize() 
    {
        Standardize1::eval(X, y, w, isd, intr, ju, xm, xs, ym, ys, xv);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        standard1_(&no, &ni, X.data(), y.data(),
                  w.data(), &i_isd, &i_intr, ju.data(),
                  xm.data(), xs.data(), &ym,
                  &ys, xv.data(), &jerr);
    }
};

struct StandardizePack : Standardize1Pack
{
    using base_t = Standardize1Pack;
    Eigen::VectorXd g;

    StandardizePack(const Eigen::MatrixXd& _X,
                    const Eigen::VectorXd& _y,
                    const Eigen::VectorXd& _w,
                    const Eigen::VectorXi& _ju,
                    bool _isd, bool _intr)
        : base_t(_X, _y, _w, _ju, _isd, _intr)
        , g(_ju.size())
    {
        g.setZero();
    }

    void standardize() 
    {
        Standardize::eval(X, y, w, isd, intr, ju, g, xm, xs, ym, ys, xv);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        standard_(&no, &ni, X.data(), y.data(),
                  w.data(), &i_isd, &i_intr, ju.data(),
                  g.data(), xm.data(), xs.data(), &ym,
                  &ys, xv.data(), &jerr);
    }
};

struct standardize_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p;
        std::tie(seed, n, p, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X(n, p);
        auto beta = dgen.make_beta(p);
        y = dgen.make_y(X, beta);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::MatrixXd X;
    Eigen::VectorXd y, w;
    Eigen::VectorXi ju;
    bool isd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_standardize1_pack(const Standardize1Pack& actual,
                                 const Standardize1Pack& expected)
    {
        expect_near_mat(actual.X, expected.X, tol);
        expect_near_vec(actual.y, expected.y, tol);
        expect_near_vec(actual.w, expected.w, tol);
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
        expect_near_vec(actual.xv, expected.xv, tol);
        EXPECT_NEAR(actual.ym, expected.ym, tol);
        EXPECT_NEAR(actual.ys, expected.ys, tol);
    }

    void check_standardize_pack(const StandardizePack& actual,
                                const StandardizePack& expected)
    {
        check_standardize1_pack(actual, expected);
        expect_near_vec(actual.g, expected.g, tol);
    }
};

TEST_P(standardize_fixture, standardize1_elnet)
{
    Standardize1Pack actual(
            X, y, w, ju, isd, intr);
    Standardize1Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_standardize1_pack(actual, expected);
}

TEST_P(standardize_fixture, standardize_elnet)
{
    StandardizePack actual(
            X, y, w, ju, isd, intr);
    StandardizePack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_standardize_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    StandardizeSuite, standardize_fixture,

    // combination of inputs: (seed, n, p, isd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool(),
        testing::Bool()
        ),
    
    // more informative test name
    [](const testing::TestParamInfo<standardize_fixture::ParamType>& info) {
      std::string name = 
          std::string("seed_") + std::to_string(std::get<0>(info.param)) + '_'
          + "n_" + std::to_string(std::get<1>(info.param)) + '_'
          + "p_" + std::to_string(std::get<2>(info.param)) + '_'
          + "isd_" + std::to_string(std::get<3>(info.param)) + '_'
          + "intr_" + std::to_string(std::get<4>(info.param));
      return name;
    }
);

// ======================================================
// Mult-response Standardize
// ======================================================

struct MultStandardize1Pack
{
    Eigen::MatrixXd X, y;
    Eigen::VectorXd w, xm, xs, ym, ys, xv;
    Eigen::VectorXi ju;
    bool isd, jsd, intr;
    double ys0 = 0;
    int jerr = 0;

    MultStandardize1Pack(
            const Eigen::MatrixXd& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _jsd, bool _intr)
        : X(_X), y(_y), w(_w), xm(_ju.size())
        , xs(_ju.size()), ym(_y.cols()), ys(_y.cols())
        , xv(_ju.size()), ju(_ju), isd(_isd), jsd(_jsd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
        ym.setZero();
        ys.setZero();
        xv.setZero();
    }

    void standardize() 
    {
        MultStandardize1::eval(X, y, w, isd, jsd, intr, ju, xm, xs, ym, ys, xv, ys0);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_jsd = jsd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        int nr = y.cols();
        multstandard1_(&no, &ni, &nr, X.data(), y.data(),
                  w.data(), &i_isd, &i_jsd, &i_intr, ju.data(),
                  xm.data(), xs.data(), ym.data(),
                  ys.data(), xv.data(), &ys0, &jerr);
    }
};

struct multstandardize1_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, bool, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p, nr;
        std::tie(seed, n, p, nr, isd, jsd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X(n, p);
        y.resize(n, nr);
        for (int i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::MatrixXd X, y;
    Eigen::VectorXd w;
    Eigen::VectorXi ju;
    bool isd, jsd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_pack(const MultStandardize1Pack& actual,
                    const MultStandardize1Pack& expected)
    {
        expect_near_mat(actual.X, expected.X, tol);
        expect_near_mat(actual.y, expected.y, tol);
        expect_near_vec(actual.w, expected.w, tol);
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
        expect_near_vec(actual.xv, expected.xv, tol);
        expect_near_vec(actual.ym, expected.ym, tol);
        expect_near_vec(actual.ys, expected.ys, tol);
        EXPECT_FLOAT_EQ(actual.ys0, expected.ys0);
    }
};

TEST_P(multstandardize1_fixture, multstandardize1_elnet)
{
    MultStandardize1Pack actual(
            X, y, w, ju, isd, jsd, intr);
    MultStandardize1Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultStandardize1Suite, multstandardize1_fixture,

    // combination of inputs: (seed, n, p, nr, isd, jsd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Values(1, 2, 3),
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

// ======================================================
// SpElnet Standardize
// ======================================================

struct SpStandardize1Pack : Standardize1Pack
{
    using base_t = Standardize1Pack;
    const Eigen::SparseMatrix<double> X;

    SpStandardize1Pack(
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _intr)
        : base_t({}, _y, _w, _ju, _isd, _intr)
        , X(_X)
    {}

    void standardize() 
    {
        SpStandardize1::eval(X, y, w, isd, intr, ju, xm, xs, ym, ys, xv);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        spstandard1_(
                &no, &ni, 
                const_cast<double*>(X.valuePtr()), 
                x_outer.data(), x_inner.data(),
                y.data(), w.data(), ju.data(), &i_isd, &i_intr,
                xm.data(), xs.data(), &ym,
                &ys, xv.data(), &jerr);
    }
};

struct SpStandardizePack : StandardizePack
{
    using base_t = StandardizePack;
    const Eigen::SparseMatrix<double> X;

    SpStandardizePack(
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _intr)
        : base_t({}, _y, _w, _ju, _isd, _intr)
        , X(_X)
    {}

    void standardize() 
    {
        SpStandardize::eval(X, y, w, isd, intr, ju, g, xm, xs, ym, ys, xv);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        spstandard_(
                &no, &ni, 
                const_cast<double*>(X.valuePtr()), 
                x_outer.data(), x_inner.data(),
                y.data(), w.data(), ju.data(), &i_isd, &i_intr,
                g.data(), xm.data(), xs.data(), &ym,
                &ys, xv.data(), &jerr);
    }
};

struct sp_standardize_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p;
        std::tie(seed, n, p, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X_sparse(n, p);
        auto beta = dgen.make_beta(p);
        y = dgen.make_y(X, beta);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::SparseMatrix<double> X;
    Eigen::VectorXd y, w;
    Eigen::VectorXi ju;
    bool isd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_standardize1_pack(const Standardize1Pack& actual,
                                 const Standardize1Pack& expected)
    {
        expect_near_vec(actual.y, expected.y, tol);
        expect_near_vec(actual.w, expected.w, tol);
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
        expect_near_vec(actual.xv, expected.xv, tol);
        EXPECT_NEAR(actual.ym, expected.ym, tol);
        EXPECT_NEAR(actual.ys, expected.ys, tol);
    }

    void check_standardize_pack(const StandardizePack& actual,
                                const StandardizePack& expected)
    {
        check_standardize1_pack(actual, expected);
        expect_near_vec(actual.g, expected.g, tol);
    }
};

TEST_P(sp_standardize_fixture, sp_standardize1_elnet)
{
    SpStandardize1Pack actual(
            X, y, w, ju, isd, intr);
    SpStandardize1Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_standardize1_pack(actual, expected);
}

TEST_P(sp_standardize_fixture, sp_standardize_elnet)
{
    SpStandardizePack actual(
            X, y, w, ju, isd, intr);
    SpStandardizePack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_standardize_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpStandardizeSuite, sp_standardize_fixture,

    // combination of inputs: (seed, n, p, isd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool(),
        testing::Bool()
        ),
    
    // more informative test name
    [](const testing::TestParamInfo<sp_standardize_fixture::ParamType>& info) {
      std::string name = 
          std::string("seed_") + std::to_string(std::get<0>(info.param)) + '_'
          + "n_" + std::to_string(std::get<1>(info.param)) + '_'
          + "p_" + std::to_string(std::get<2>(info.param)) + '_'
          + "isd_" + std::to_string(std::get<3>(info.param)) + '_'
          + "intr_" + std::to_string(std::get<4>(info.param));
      return name;
    }
);

// ======================================================
// Mult-response SpStandardize
// ======================================================

struct MultSpStandardize1Pack
{
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd y;
    Eigen::VectorXd w, xm, xs, ym, ys, xv;
    Eigen::VectorXi ju;
    bool isd, jsd, intr;
    double ys0 = 0;
    int jerr = 0;

    MultSpStandardize1Pack(
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::MatrixXd& _y,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _jsd, bool _intr)
        : X(_X), y(_y), w(_w), xm(_ju.size())
        , xs(_ju.size()), ym(_y.cols()), ys(_y.cols())
        , xv(_ju.size()), ju(_ju), isd(_isd), jsd(_jsd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
        ym.setZero();
        ys.setZero();
        xv.setZero();
    }

    void standardize() 
    {
        MultSpStandardize1::eval(X, y, w, isd, jsd, intr, ju, xm, xs, ym, ys, xv, ys0);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_jsd = jsd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        int nr = y.cols();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        multspstandard1_(&no, &ni, &nr, 
                X.valuePtr(), x_outer.data(), x_inner.data(), y.data(),
                w.data(), ju.data(), &i_isd, &i_jsd, &i_intr, 
                xm.data(), xs.data(), ym.data(),
                ys.data(), xv.data(), &ys0, &jerr);
    }
};

struct multspstandardize1_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, bool, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p, nr;
        std::tie(seed, n, p, nr, isd, jsd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X_sparse(n, p);
        y.resize(n, nr);
        for (int i = 0; i < nr; ++i) {
            auto beta = dgen.make_beta(p);
            y.col(i) = dgen.make_y(X, beta);
        }
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::SparseMatrix<double> X;
    Eigen::MatrixXd y;
    Eigen::VectorXd w;
    Eigen::VectorXi ju;
    bool isd, jsd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_pack(const MultSpStandardize1Pack& actual,
                    const MultSpStandardize1Pack& expected)
    {
        expect_near_mat(actual.y, expected.y, tol);
        expect_near_vec(actual.w, expected.w, tol);
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
        expect_near_vec(actual.xv, expected.xv, tol);
        expect_near_vec(actual.ym, expected.ym, tol);
        expect_near_vec(actual.ys, expected.ys, tol);
        EXPECT_FLOAT_EQ(actual.ys0, expected.ys0);
    }
};

TEST_P(multspstandardize1_fixture, multspstandardize1_elnet)
{
    MultSpStandardize1Pack actual(
            X, y, w, ju, isd, jsd, intr);
    MultSpStandardize1Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultSpStandardize1Suite, multspstandardize1_fixture,

    // combination of inputs: (seed, n, p, nr, isd, jsd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Values(1, 2, 3),
        testing::Bool(),
        testing::Bool(),
        testing::Bool()
        )
);

// ======================================================
// Lognet Standardize
// ======================================================

struct LStandardize1Pack
{
    Eigen::MatrixXd X;
    const Eigen::VectorXd& w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd xm, xs;
    bool isd, intr;

    LStandardize1Pack(
            const Eigen::MatrixXd& _X,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _intr)
        : X(_X), w(_w), ju(_ju), xm(_ju.size())
        , xs(_ju.size()), isd(_isd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
    }

    void standardize() 
    {
        LStandardize1::eval(X, w, ju, isd, intr, xm, xs);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        lstandard1_(&no, &ni, X.data(),
                  const_cast<double*>(w.data()), 
                  const_cast<int*>(ju.data()), &i_isd, &i_intr,
                  xm.data(), xs.data());
    }
};

struct lstandardize1_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p;
        std::tie(seed, n, p, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X(n, p);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::MatrixXd X;
    Eigen::VectorXd w;
    Eigen::VectorXi ju;
    bool isd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_standardize1_pack(const LStandardize1Pack& actual,
                                 const LStandardize1Pack& expected)
    {
        expect_near_mat(actual.X, expected.X, tol);
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
    }
};

TEST_P(lstandardize1_fixture, lstandardize1_lognet)
{
    LStandardize1Pack actual(
            X, w, ju, isd, intr);
    LStandardize1Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_standardize1_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    LStandardize1Suite, lstandardize1_fixture,

    // combination of inputs: (seed, n, p, isd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool(),
        testing::Bool()
        ),
    
    // more informative test name
    [](const testing::TestParamInfo<standardize_fixture::ParamType>& info) {
      std::string name = 
          std::string("seed_") + std::to_string(std::get<0>(info.param)) + '_'
          + "n_" + std::to_string(std::get<1>(info.param)) + '_'
          + "p_" + std::to_string(std::get<2>(info.param)) + '_'
          + "isd_" + std::to_string(std::get<3>(info.param)) + '_'
          + "intr_" + std::to_string(std::get<4>(info.param));
      return name;
    }
);

// ======================================================
// Multi-class Group Lasso Standardize
// ======================================================

struct MultLStandardize1Pack
{
    Eigen::MatrixXd X;
    const Eigen::VectorXd& w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd xm, xs, xv;
    bool isd, intr;

    MultLStandardize1Pack(
            const Eigen::MatrixXd& _X,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _intr)
        : X(_X), w(_w), ju(_ju), xm(_ju.size())
        , xs(_ju.size()), xv(_ju.size()), isd(_isd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
        xv.setZero();
    }

    void standardize() 
    {
        MultLStandardize1::eval(X, w, ju, isd, intr, xm, xs, xv);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        multlstandard1_(&no, &ni, X.data(),
                  const_cast<double*>(w.data()), 
                  const_cast<int*>(ju.data()), &i_isd, &i_intr,
                  xm.data(), xs.data(), xv.data());
    }
};

struct multlstandardize1_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p;
        std::tie(seed, n, p, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X(n, p);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::MatrixXd X;
    Eigen::VectorXd w;
    Eigen::VectorXi ju;
    bool isd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_pack(const MultLStandardize1Pack& actual,
                    const MultLStandardize1Pack& expected)
    {
        expect_near_mat(actual.X, expected.X, tol);
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
        expect_near_vec(actual.xv, expected.xv, tol);
    }
};

TEST_P(multlstandardize1_fixture, multlstandardize1_lognet)
{
    MultLStandardize1Pack actual(
            X, w, ju, isd, intr);
    MultLStandardize1Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultLStandardize1Suite, multlstandardize1_fixture,

    // combination of inputs: (seed, n, p, isd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool(),
        testing::Bool()
        )
);

// ======================================================
// SpLognet Standardize
// ======================================================

struct SpLStandardize2Pack
{
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd xm, xs;
    bool isd, intr;

    SpLStandardize2Pack(
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _intr)
        : X(_X), w(_w), ju(_ju), xm(_ju.size())
        , xs(_ju.size()), isd(_isd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
    }

    void standardize() 
    {
        SpLStandardize2::eval(X, w, ju, isd, intr, xm, xs);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        splstandard2_(&no, &ni, 
                const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(),
                const_cast<double*>(w.data()), 
                const_cast<int*>(ju.data()), &i_isd, &i_intr,
                xm.data(), xs.data());
    }
};

struct sp_lstandardize2_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p;
        std::tie(seed, n, p, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X_sparse(n, p);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::SparseMatrix<double> X;
    Eigen::VectorXd w;
    Eigen::VectorXi ju;
    bool isd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_pack(const SpLStandardize2Pack& actual,
                    const SpLStandardize2Pack& expected)
    {
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
    }
};

TEST_P(sp_lstandardize2_fixture, sp_lstandardize2_lognet)
{
    SpLStandardize2Pack actual(
            X, w, ju, isd, intr);
    SpLStandardize2Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpLStandardize2Suite, sp_lstandardize2_fixture,

    // combination of inputs: (seed, n, p, isd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool(),
        testing::Bool()
        )
);

// ======================================================
// MultSpLognet Standardize
// ======================================================

struct MultSpLStandardize2Pack
{
    const Eigen::SparseMatrix<double>& X;
    const Eigen::VectorXd& w;
    const Eigen::VectorXi& ju;
    Eigen::VectorXd xm, xs, xv;
    bool isd, intr;

    MultSpLStandardize2Pack(
            const Eigen::SparseMatrix<double>& _X,
            const Eigen::VectorXd& _w,
            const Eigen::VectorXi& _ju,
            bool _isd, bool _intr)
        : X(_X), w(_w), ju(_ju), xm(_ju.size())
        , xs(_ju.size()), xv(_ju.size()), isd(_isd), intr(_intr)
    {
        xm.setZero();
        xs.setZero();
    }

    void standardize() 
    {
        MultSpLStandardize2::eval(X, w, ju, isd, intr, xm, xs, xv);
    }

    void standardize_legacy()
    {
        int i_isd = isd;
        int i_intr = intr;
        int no = X.rows();
        int ni = X.cols();
        auto x_inner = make_sp_inner_idx_1idx(X);
        auto x_outer = make_sp_outer_idx_1idx(X);
        multsplstandard2_(&no, &ni, 
                const_cast<double*>(X.valuePtr()), x_outer.data(), x_inner.data(),
                const_cast<double*>(w.data()), 
                const_cast<int*>(ju.data()), &i_isd, &i_intr,
                xm.data(), xs.data(), xv.data());
    }
};

struct multsp_lstandardize2_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool, bool>
    >
{
    void SetUp() override
    {
        int seed, n, p;
        std::tie(seed, n, p, isd, intr) = GetParam();
        DataGen dgen(seed);
        X = dgen.make_X_sparse(n, p);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
    }

protected: 
    Eigen::SparseMatrix<double> X;
    Eigen::VectorXd w;
    Eigen::VectorXi ju;
    bool isd, intr;

    // double precision error: 
    // absolute difference shouldn't be any more than this amount
    static constexpr double tol = 1e-14;

    void check_pack(const MultSpLStandardize2Pack& actual,
                    const MultSpLStandardize2Pack& expected)
    {
        expect_near_vec(actual.xm, expected.xm, tol);
        expect_near_vec(actual.xs, expected.xs, tol);
        expect_near_vec(actual.xv, expected.xv, tol);
    }
};

TEST_P(multsp_lstandardize2_fixture, multsp_lstandardize2_lognet)
{
    MultSpLStandardize2Pack actual(
            X, w, ju, isd, intr);
    MultSpLStandardize2Pack expected(actual);
    actual.standardize();
    expected.standardize_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    MultSpLStandardize2Suite, multsp_lstandardize2_fixture,

    // combination of inputs: (seed, n, p, isd, intr)
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool(),
        testing::Bool()
        )
);

} // namespace glmnetpp
