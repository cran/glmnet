#include "gtest/gtest.h"
#include <testutil/data_util.hpp>
#include <legacy/legacy.h>
#include <glmnetpp_bits/internal.hpp>

namespace glmnetpp {

struct InternalPack
{
    double sml, eps, big, rsqmax, pmin, exmx;
    int mnlam, itrace;

    void get_int_parms()
    {
        glmnetpp::get_int_parms(sml, eps, big, mnlam, rsqmax, pmin, exmx, itrace);
    }
    
    void get_int_parms_()
    {
        ::get_int_parms_(&sml, &eps, &big, &mnlam, &rsqmax, &pmin, &exmx, &itrace);
    }
};

struct internal_fixture : ::testing::Test
{
protected:

    void check_internal() const
    {
        InternalPack actual, expected;
        actual.get_int_parms();
        expected.get_int_parms_();
        EXPECT_DOUBLE_EQ(actual.sml, expected.sml);
        EXPECT_DOUBLE_EQ(actual.eps, expected.eps);
        EXPECT_DOUBLE_EQ(actual.big, expected.big);

        // due to single-precision floating-pt error, fortran returns 0.999 + O(1e-8)
        EXPECT_NEAR(actual.rsqmax, expected.rsqmax, 1e-7);

        EXPECT_DOUBLE_EQ(actual.pmin, expected.pmin);
        EXPECT_DOUBLE_EQ(actual.exmx, expected.exmx);
        EXPECT_EQ(actual.mnlam, expected.mnlam);
        EXPECT_EQ(actual.itrace, expected.itrace);
    }
};

TEST_F(internal_fixture, get_int_parms_compat)
{
    check_internal();
}

template <class ArgType>
struct chg_fixture : 
    internal_fixture,
    ::testing::WithParamInterface<ArgType>
{
protected:
};

using chg_fixture_double = chg_fixture<double>;
using chg_fixture_int = chg_fixture<int>;

#define GENERATE_CHG_TEST(fixture, chg_name) \
    TEST_P(fixture, chg_name##_compat) \
    { \
        auto arg = GetParam(); \
        chg_name(arg); \
        chg_name##_(&arg); \
        check_internal();  \
    }

GENERATE_CHG_TEST(chg_fixture_double, chg_fract_dev)
GENERATE_CHG_TEST(chg_fixture_double, chg_min_flmin)
GENERATE_CHG_TEST(chg_fixture_double, chg_dev_max)
GENERATE_CHG_TEST(chg_fixture_double, chg_big)
GENERATE_CHG_TEST(chg_fixture_double, chg_min_null_prob)
GENERATE_CHG_TEST(chg_fixture_double, chg_max_exp)
GENERATE_CHG_TEST(chg_fixture_int, chg_min_lambdas)
GENERATE_CHG_TEST(chg_fixture_int, chg_itrace)

INSTANTIATE_TEST_SUITE_P(
    ChgDoubleSuite, chg_fixture_double,
    testing::Values(1., 242., 1.52, -92.2, 15.2, 23.4)
);

INSTANTIATE_TEST_SUITE_P(
    ChgIntSuite, chg_fixture_int,
    testing::Values(-1, 0, 5, 2, 199, 203284, 23, -238, -52, 32)
);

} // namespace glmnetpp
