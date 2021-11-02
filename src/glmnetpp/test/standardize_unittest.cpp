#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/standardize.hpp>
#include <tuple>

namespace glmnetpp {

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

// =================================================================================
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

} // namespace glmnetpp
