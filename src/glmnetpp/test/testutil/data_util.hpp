#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace glmnetpp {

inline Eigen::MatrixXd read_csv(const std::string& filename, char delim=',')
{
    std::vector<double> matrixEntries;

    std::ifstream matrixDataFile(filename);

    std::string matrixRowString;
    std::string matrixEntry;

    int matrixRowNumber = 0;

    while (std::getline(matrixDataFile, matrixRowString)) 
    {
        std::stringstream matrixRowStringStream(matrixRowString); 
        while (std::getline(matrixRowStringStream, matrixEntry, delim)) 
        {
            matrixEntries.push_back(stod(matrixEntry));
        }
        matrixRowNumber++;
    }

    // here we convet the vector variable into the matrix and return the resulting object,
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}

// define csv format
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

template <typename Derived>
void write_csv(const std::string& name, const Eigen::MatrixBase<Derived>& matrix)
{
    std::ofstream file(name.c_str());
    file << std::setprecision(16) << matrix.format(CSVFormat);
}

// Center and standardize a matrix with 0 mean columns and scaled by var/n each column
template <class T, int R, int C>
inline Eigen::Matrix<T, R, C> 
center_scale(const Eigen::Matrix<T, R, C>& X)
{
    Eigen::Matrix<T, R, C> out(X.rows(), X.cols());
    auto n = X.rows();
    for (int i = 0; i < X.cols(); ++i) {
        out.col(i) = X.col(i).array() - X.col(i).mean();
        out.col(i) /= out.col(i).norm() / std::sqrt(n);
    }
    return out;
}

struct DataGen
{

    DataGen(size_t seed)
        : gen(seed)
    {}

    auto make_X(size_t n, size_t p) 
    {
        std::normal_distribution<double> norm(0., 1.);
        Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(
                n, p, [&](auto, auto) { return norm(gen); });
        return X;
    }

    auto make_X_sparse(size_t n, size_t p, double sparsity=0.4)
    {
        Eigen::SparseMatrix<double> X(n, p);
        std::normal_distribution<double> norm(0., 1.);
        std::bernoulli_distribution bern(sparsity);
        for (size_t j = 0; j < p; ++j) {
            X.coeffRef(0, j) = norm(gen);   // always make sure first row is non-zero so that stddev > 0
            for (size_t i = 1; i < n; ++i) {
                if (bern(gen)) X.coeffRef(i, j) = norm(gen);
            }
        }
        X.makeCompressed();
        return X;
    }

    auto make_beta(size_t p, double sparsity=0.5) 
    {
        std::bernoulli_distribution bern_sp(sparsity);
        std::normal_distribution<double> norm(0., 1.);
        Eigen::VectorXd beta = Eigen::VectorXd::NullaryExpr(
                p, [&](auto) { return bern_sp(gen) * norm(gen); });
        return beta;
    }

    auto make_y(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) 
    {
        std::normal_distribution<double> norm(0., 1.);
        Eigen::VectorXd y;
        auto n = X.rows();
        y = X * beta + y.NullaryExpr(n, [&](auto) { return norm(gen); });
        auto ym = y.mean();
        auto ys = std::sqrt((y.array() - y.mean()).square().sum());
        y.array() -= ym;
        y /= ys;
        return y;
    }

    auto make_w(size_t n) 
    {
        std::bernoulli_distribution bern_half(0.5);
        Eigen::VectorXd weights;
        weights = weights.NullaryExpr(n, [&](auto) { return bern_half(gen) + 1; });
        weights /= weights.sum();
        return weights;
    }

    auto make_ju(size_t p, double inclusion_rate=0.99) 
    {
        std::bernoulli_distribution bern_in(inclusion_rate);
        Eigen::VectorXi inclusion;
        inclusion = inclusion.NullaryExpr(p, [&](auto) { return bern_in(gen); });
        return inclusion;
    }

    auto make_jd(size_t p)
    {
        std::uniform_int_distribution<> mult(1, p);
        auto n_jd = mult(gen);
        Eigen::VectorXi jd(n_jd+1);
        jd = Eigen::VectorXi::NullaryExpr(n_jd+1, [&](auto){ return mult(gen); });
        jd[0] = n_jd;
        return jd;
    }

    auto make_vp(size_t p)   
    {
        std::uniform_int_distribution<> mult(1, p);
        Eigen::VectorXd vp;
        vp = vp.NullaryExpr(p, [&](auto) { return mult(gen); });
        vp /= vp.sum() / p;
        return vp;
    }

    auto make_cl(size_t p) 
    {
        std::normal_distribution<double> norm(0., 1.);
        Eigen::MatrixXd cl;
        cl = Eigen::MatrixXd::NullaryExpr(2, p, 
            [&](auto i, auto) {
                if (i == 0) return -100 * std::abs(norm(gen));
                else return 100 * std::abs(norm(gen));
            });
        return cl;
    }

    auto make_nx(size_t p) 
    {
        std::uniform_int_distribution<> mult(1, p);
        return mult(gen);
    }

    auto make_ne(size_t p)
    {
        std::uniform_int_distribution<> mult(1, p);
        return mult(gen);
    }

    auto make_ulam(size_t nlam) 
    {
        std::uniform_real_distribution<double> unif(0.0001, 2.);
        Eigen::VectorXd ulam;
        ulam = ulam.NullaryExpr(nlam, [&](auto){ return unif(gen); });
        std::sort(ulam.data(), ulam.data() + ulam.size(), std::greater<double>());
        return ulam;
    }

private:
    std::mt19937 gen;
};


// Make sparse matrix inner index array 1-indexed.
template <class SparseMatType>
inline auto make_sp_inner_idx_1idx(const SparseMatType& X)
{
    Eigen::VectorXi x_inner(X.nonZeros());
    for (int i = 0; i < x_inner.size(); ++i) {
        x_inner[i] = X.innerIndexPtr()[i] + 1;
    }
    return x_inner;
}

// Make sparse matrix outer index array 1-indexed.
template <class SparseMatType>
inline auto make_sp_outer_idx_1idx(const SparseMatType& X)
{
    Eigen::VectorXi x_outer(X.cols() + 1);
    for (int i = 0; i < x_outer.size(); ++i) {
        x_outer[i] = X.outerIndexPtr()[i] + 1;
    }
    return x_outer;
}

} // namespace glmnetpp
