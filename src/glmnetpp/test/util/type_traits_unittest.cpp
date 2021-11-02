#include "gtest/gtest.h"
#include <glmnetpp_bits/util/type_traits.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace glmnetpp {
namespace util {

// Pair is expected to be a tuple<T, bool>
// where T is the type to test and bool is true if T is expected to be dense.
template <class Pair>
struct is_dense_fixture : ::testing::Test
{};

using list_t = ::testing::Types<
    std::tuple<double, std::false_type>,
    std::tuple<float, std::false_type>,
    std::tuple<int, std::false_type>,
    std::tuple<Eigen::SparseMatrix<double>, std::false_type>,
    std::tuple<Eigen::Map<Eigen::SparseMatrix<double>>, std::false_type>,
    std::tuple<Eigen::Map<Eigen::SparseMatrix<int>>, std::false_type>,
    std::tuple<Eigen::MappedSparseMatrix<double>, std::false_type>,
    std::tuple<Eigen::VectorXd, std::true_type>,
    std::tuple<Eigen::MatrixXd, std::true_type>,
    std::tuple<Eigen::VectorXi, std::true_type>,
    std::tuple<Eigen::MatrixXi, std::true_type>,
    std::tuple<Eigen::ArrayXd, std::true_type>,
    std::tuple<Eigen::ArrayXi, std::true_type>,
    std::tuple<Eigen::Matrix3f, std::true_type>,
    std::tuple<Eigen::Matrix3d, std::true_type>
>;

TYPED_TEST_SUITE(is_dense_fixture, list_t,);

template <bool dense, class T>
struct CheckDense
{
    static_assert(is_dense<T>::value, "Dense check failed.");
};

template <class T>
struct CheckDense<false, T>
{
    static_assert(!is_dense<T>::value, "Sparse check failed.");
};

TYPED_TEST(is_dense_fixture, 
            test_double)
{
    using T = std::tuple_element_t<0, TypeParam>;
    using is_dense_type = std::tuple_element_t<1, TypeParam>;
    CheckDense<std::is_same<is_dense_type, std::true_type>::value, T> dummy;
    static_cast<void>(dummy);
}

} // namespace util
} // namespace glmnetpp
