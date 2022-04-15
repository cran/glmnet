#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/elnet_driver/chkvars.hpp>
#include <Eigen/SparseCore>
#include <random>
#include <tuple>

namespace glmnetpp {

struct ChkvarsPack
{
    const Eigen::MatrixXd& X;
    Eigen::VectorXi ju;
    ChkvarsPack(const Eigen::MatrixXd& _X)
        : X(_X), ju(X.cols())
    {
        ju.setZero();
    }

    void chkvars()
    {
        Chkvars::eval(X, ju);
    }

    void chkvars_legacy()
    {
        int no = X.rows();
        int ni = X.cols();
        chkvars_(&no, &ni, const_cast<double*>(X.data()), ju.data());
    }
};

struct chkvars_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool>
    >
{
protected: 

    auto gen_data(size_t seed, size_t n, size_t p, bool force_same_value)
    {
        std::mt19937 gen;
        gen.seed(seed);
        std::normal_distribution<double> norm(0.,1.);
        Eigen::MatrixXd X;
        X = X.NullaryExpr(n, p, [&](auto, auto) { return norm(gen); });
        if (force_same_value) {
            std::bernoulli_distribution bern_half(0.5);
            for (size_t j = 0; j < p; ++j) {
                if (bern_half(gen)) {
                    X.col(j).array() = X(0, j);
                }
            }
        }
        return X;
    }

    void check_pack(const ChkvarsPack& actual,
                    const ChkvarsPack& expected)
    {
        expect_eq_vec(actual.ju, expected.ju);
    }
};

TEST_P(chkvars_fixture, chkvars_elnet)
{
    int seed, n, p;
    bool force_same_value;
    std::tie(seed, n, p, force_same_value) = GetParam();
    auto X = gen_data(seed, n, p, force_same_value);
    ChkvarsPack actual(X), expected(X);
    actual.chkvars();
    expected.chkvars_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    ChkvarsSuite, chkvars_fixture,
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool()
    ),

    [](const testing::TestParamInfo<chkvars_fixture::ParamType>& info) {
      std::string name = 
          std::string("seed_") + std::to_string(std::get<0>(info.param)) + '_'
          + "n_" + std::to_string(std::get<1>(info.param)) + '_'
          + "p_" + std::to_string(std::get<2>(info.param)) + '_'
          + "force_" + std::to_string(std::get<3>(info.param));
      return name;
    }
);

// ========================================================

struct SpChkvarsPack
{
    const Eigen::SparseMatrix<double>& X;
    Eigen::VectorXi ju;

    SpChkvarsPack(const Eigen::SparseMatrix<double>& _X)
        : X(_X), ju(X.cols())
    {
        ju.setZero();
    }

    void chkvars()
    {
        SpChkvars::eval(X, ju);
    }

    void chkvars_legacy()
    {
        int no = X.rows();
        int ni = X.cols();
        spchkvars_(&no, &ni, 
                   const_cast<double*>(X.valuePtr()), 
                   const_cast<int*>(X.outerIndexPtr()),
                   ju.data());
    }
};

struct sp_chkvars_fixture : 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, bool>
    >
{
protected: 

    auto gen_data(size_t seed, size_t n, size_t p, bool force_same_value)
    {
        DataGen dgen(seed);
        auto X = dgen.make_X_sparse(n, p);
        if (force_same_value) {
            std::mt19937 gen(seed);
            std::bernoulli_distribution bern_half(0.5);
            for (int j = 0; j < X.cols(); ++j) {
                if (bern_half(gen)) {
                    auto t = X.coeff(0,j);
                    for (int i = 0; i < X.rows(); ++i) {
                        X.coeffRef(i, j) = t;
                    }
                }
            }
        }
        X.makeCompressed();
        return X;
    }

    void check_pack(const SpChkvarsPack& actual,
                    const SpChkvarsPack& expected)
    {
        expect_eq_vec(actual.ju, expected.ju);
    }
};

TEST_P(sp_chkvars_fixture, sp_chkvars_elnet)
{
    int seed, n, p;
    bool force_same_value;
    std::tie(seed, n, p, force_same_value) = GetParam();
    auto X = gen_data(seed, n, p, force_same_value);
    SpChkvarsPack actual(X), expected(X);
    actual.chkvars();
    expected.chkvars_legacy();
    check_pack(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(
    SpChkvarsSuite, sp_chkvars_fixture,
    testing::Combine(
        testing::Values(10, 23, 145, 241, 412, 23968, 31),
        testing::Values(10, 20, 30, 50),
        testing::Values(10, 20, 30, 40),
        testing::Bool()
    ),

    [](const testing::TestParamInfo<chkvars_fixture::ParamType>& info) {
      std::string name = 
          std::string("seed_") + std::to_string(std::get<0>(info.param)) + '_'
          + "n_" + std::to_string(std::get<1>(info.param)) + '_'
          + "p_" + std::to_string(std::get<2>(info.param)) + '_'
          + "force_" + std::to_string(std::get<3>(info.param));
      return name;
    }
);

} // namespace glmnetpp
