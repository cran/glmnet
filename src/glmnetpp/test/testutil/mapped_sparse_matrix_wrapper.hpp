#pragma once
#include <Eigen/Core>

namespace glmnetpp {

// A wrapper class to represent a standardized sparse (column) vector.
template <class ColVecType>
struct MappedSparseVectorStandardized
    : public ColVecType
{
private:
    using base_t = ColVecType;
    using value_t = typename base_t::Scalar; 
    
public:    
    MappedSparseVectorStandardized(
        const base_t& c,
        value_t mean,
        value_t sd)
        : base_t(c)
        , mean_(mean)
        , sd_(sd)
    {}
    
    value_t mean() const { return mean_; }
    value_t sd() const { return sd_; }
    
private:
    const value_t mean_;
    const value_t sd_;
};

// A wrapper class to represent a standardized sparse matrix.
template <class SparseMatMapType>
struct MappedSparseMatrixWrapper
    : public SparseMatMapType
{   
private:
    using base_t = SparseMatMapType;
    using value_t = typename SparseMatMapType::Scalar; 
    using vec_t = Eigen::Matrix<value_t, Eigen::Dynamic, 1>;
    using map_vec_t = Eigen::Map<vec_t>;
    using index_t = typename SparseMatMapType::Index;
                                
public:
    MappedSparseMatrixWrapper(
        const base_t& m, 
        const map_vec_t& mean,
        const map_vec_t& sd)
        : base_t(m)
        , mean_(mean)
        , sd_(sd)
    {} 
    
    auto col(index_t j) const { 
        auto col_j = base_t::col(j);
        using col_t = std::decay_t<decltype(col_j)>;
        return MappedSparseVectorStandardized<col_t>(
            col_j, mean_(j), sd_(j));
    }
    
    auto col(index_t j) { 
        auto col_j = base_t::col(j);
        using col_t = std::decay_t<decltype(col_j)>;
        return MappedSparseVectorStandardized<col_t>(
            col_j, mean_(j), sd_(j));
    }
    
private:
    const map_vec_t mean_;
    const map_vec_t sd_;
};

template <class SparseMatMapType
        , class MeanType
        , class SDType>
inline auto 
make_mapped_sparse_matrix_wrapper(
        const SparseMatMapType& m,
        const MeanType& mean,
        const SDType& sd)
{ return MappedSparseMatrixWrapper<SparseMatMapType>(m, mean, sd); }

} // namespace glmnetpp
