#pragma once
#include "gtest/gtest.h"
#include <Eigen/Core>

namespace glmnetpp {

struct base_fixture : ::testing::Test
{
protected:
    using value_t = double;
    using index_t = Eigen::Index;

    // Useful tools to test vector equality
#define expect_double_eq_vec(v1, v2) \
    { \
        EXPECT_EQ(v1.size(), v2.size()); \
        for (index_t i = 0; i < v1.size(); ++i) { \
            EXPECT_DOUBLE_EQ(v1[i], v2[i]); \
        } \
    } 

#define expect_float_eq_vec(v1, v2) \
    { \
        EXPECT_EQ(v1.size(), v2.size()); \
        for (index_t i = 0; i < v1.size(); ++i) { \
            EXPECT_FLOAT_EQ( \
                    static_cast<float>(v1[i]), \
                    static_cast<float>(v2[i]) \
                    ); \
        } \
    } 

#define expect_eq_vec(v1, v2) \
    { \
        EXPECT_EQ(v1.size(), v2.size()); \
        for (index_t i = 0; i < v1.size(); ++i) { \
            EXPECT_EQ(v1[i], v2[i]); \
        } \
    }

#define expect_double_eq_mat(m1, m2) \
    { \
        EXPECT_EQ(m1.rows(), m2.rows()); \
        EXPECT_EQ(m1.cols(), m2.cols()); \
        for (index_t j = 0; j < m1.cols(); ++j) { \
            for (index_t i = 0; i < m1.rows(); ++i) { \
                EXPECT_DOUBLE_EQ(m1(i,j), m2(i,j)); \
            } \
        } \
    }

#define expect_float_eq_mat(m1, m2) \
    { \
        EXPECT_EQ(m1.rows(), m2.rows()); \
        EXPECT_EQ(m1.cols(), m2.cols()); \
        for (index_t j = 0; j < m1.cols(); ++j) { \
            for (index_t i = 0; i < m1.rows(); ++i) { \
                EXPECT_FLOAT_EQ( \
                    static_cast<float>(m1(i,j)), \
                    static_cast<float>(m2(i,j))); \
            } \
        } \
    }
    
#define expect_near_vec(v1, v2, tol) \
    { \
        EXPECT_EQ(v1.size(), v2.size()); \
        for (index_t i = 0; i < v1.size(); ++i) { \
            EXPECT_NEAR(v1[i], v2[i], tol); \
        } \
    }

#define expect_near_mat(m1, m2, tol) \
    { \
        EXPECT_EQ(m1.rows(), m2.rows()); \
        EXPECT_EQ(m1.cols(), m2.cols()); \
        for (index_t j = 0; j < m1.cols(); ++j) { \
            for (index_t i = 0; i < m1.rows(); ++i) { \
                EXPECT_NEAR(m1(i,j), m2(i,j), tol); \
            } \
        } \
    }
};

} // namespace glmnetpp
