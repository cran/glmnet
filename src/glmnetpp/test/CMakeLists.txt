# Try to find package Threads for the use of pthreads
unset(GLMNETPP_HAS_PTHREAD)
find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
    set(GLMNETPP_HAS_PTHREAD ON)
endif()

set(GLMNETPP_HAS_PTHREAD_MACRO "")
if (DEFINED GLMNETPP_HAS_PTHREAD)
    set(GLMNETPP_HAS_PTHREAD_MACRO -DGLMNETPP_HAS_PTHREAD)
endif()

set(GLMNETPP_TEST_COMPILE_FLAGS ${GLMNETPP_STRICT_WARNINGS})
set(GLMNETPP_TEST_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${GLMNETPP_SOURCEDIR})
set(GLMNETPP_TEST_LIBS ${PROJECT_NAME} GTest::gtest_main Eigen3::Eigen)
set(GLMNETPP_TEST_COMPILE_DEFS "")

# Append variables depending on compiler
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(GLMNETPP_TEST_COMPILE_FLAGS ${GLMNETPP_TEST_COMPILE_FLAGS} -fopenmp)
    set(GLMNETPP_TEST_LIBS ${GLMNETPP_TEST_LIBS} -fopenmp)
endif()

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(GLMNETPP_TEST_LIBS ${GLMNETPP_TEST_LIBS} pthread)
endif()

if (GLMNETPP_ENABLE_COVERAGE)
    set(GLMNETPP_TEST_LIBS ${GLMNETPP_TEST_LIBS} gcov)
endif()

# This is if we can't compile Fortran code,
# but we still want to run our C++ method.
# We enable the mocked version of the legacy Fortran code.
if (GLMNETPP_MOCK_LEGACY)
    set(GLMNETPP_TEST_COMPILE_DEFS ${GLMNETPP_TEST_COMPILE_DEFS} -DGLMNETPP_MOCK_LEGACY)
endif()

########################################################################
# Utility TEST
########################################################################

add_executable(utility_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/util/type_traits_unittest.cpp
    )
target_compile_options(utility_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_FLAGS})
target_include_directories(utility_unittest PRIVATE ${GLMNETPP_TEST_INCLUDES})
target_link_libraries(utility_unittest PRIVATE ${GLMNETPP_TEST_LIBS})
add_test(utility_unittest utility_unittest)

########################################################################
# Translation TEST
########################################################################

add_executable(translation_unittest
    # gaussian
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet1_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet2_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multelnet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multelnet2_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/spelnet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/spelnet1_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/spelnet2_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multspelnet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multspelnet2_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/wls_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/wls_dscheck_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet_dscheck_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multelnet_dscheck_unittest.cpp

    # binomial
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/lognet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/lognet2n_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/lognetn_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multlognetn_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/splognet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/splognet2n_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/splognetn_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/multsplognetn_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/lognet_dscheck_unittest.cpp

    # poisson
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/fishnet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/spfishnet_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/fishnet1_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/spfishnet1_unittest.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/translation/fishnet_dscheck_unittest.cpp

    # Must link with utility files
    ${CMAKE_CURRENT_SOURCE_DIR}/testutil/internal.cpp

    # Explicitly add Fortran source for debugging purposes
    ${GLMNETPP_SOURCEDIR}/legacy/glmnet5dpclean.f
    ${GLMNETPP_SOURCEDIR}/legacy/wls.f
    ${GLMNETPP_SOURCEDIR}/legacy/pb.c   # only needed for fortran code to link to this
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # DO NOT put -Werror because Fortran gives a bunch of warnings that we unfortunately cannot change.
    target_compile_options(translation_unittest PRIVATE -Wextra -Wpedantic)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(translation_unittest PRIVATE -fopenmp)
endif()
target_compile_options(translation_unittest PRIVATE -g -Wall)
target_compile_definitions(translation_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_DEFS})
target_include_directories(translation_unittest PRIVATE
    ${GLMNETPP_TEST_INCLUDES}
    ${R_INCLUDE_DIR})
target_link_libraries(translation_unittest PRIVATE ${GLMNETPP_TEST_LIBS})
if (NOT GLMNETPP_MOCK_LEGACY)
    target_link_libraries(translation_unittest PRIVATE ${RLIB})
endif()
add_test(translation_unittest translation_unittest)

########################################################################
# Elnet Point TEST
########################################################################
add_executable(elnet_point_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_point/gaussian_wls_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_point/sp_gaussian_wls_unittest.cpp
    )
target_compile_options(elnet_point_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_FLAGS})
target_compile_definitions(elnet_point_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_DEFS})
target_include_directories(elnet_point_unittest PRIVATE ${GLMNETPP_TEST_INCLUDES})
target_link_libraries(elnet_point_unittest PRIVATE ${GLMNETPP_TEST_LIBS})
add_test(elnet_point_unittest elnet_point_unittest)

########################################################################
# Elnet Path TEST
########################################################################

add_executable(elnet_path_unittest
    # gaussian
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/gaussian_cov_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/gaussian_naive_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/gaussian_multi_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_gaussian_cov_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_gaussian_naive_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_gaussian_multi_unittest.cpp
    # binomial
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/binomial_two_class_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/binomial_multi_class_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/binomial_multi_class_group_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_binomial_two_class_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_binomial_multi_class_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_binomial_multi_class_group_unittest.cpp
    # poisson
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/poisson_naive_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_poisson_naive_unittest.cpp
    # Must link with utility files
    ${CMAKE_CURRENT_SOURCE_DIR}/testutil/internal.cpp
    )

target_compile_options(elnet_path_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_FLAGS})
target_compile_definitions(elnet_path_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_DEFS})
target_include_directories(elnet_path_unittest PRIVATE ${GLMNETPP_TEST_INCLUDES})
target_link_libraries(elnet_path_unittest PRIVATE ${GLMNETPP_TEST_LIBS})
add_test(elnet_path_unittest elnet_path_unittest)

########################################################################
# Elnet Driver TEST (Integration Test)
########################################################################

add_executable(elnet_driver_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_driver/gaussian_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_driver/binomial_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_driver/poisson_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_driver/standardize_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_driver/chkvars_unittest.cpp
    # Must link with utility files
    ${CMAKE_CURRENT_SOURCE_DIR}/testutil/internal.cpp
    )

target_compile_options(elnet_driver_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_FLAGS})
target_compile_definitions(elnet_driver_unittest PRIVATE ${GLMNETPP_TEST_COMPILE_DEFS})
target_include_directories(elnet_driver_unittest PRIVATE ${GLMNETPP_TEST_INCLUDES})
target_link_libraries(elnet_driver_unittest PRIVATE ${GLMNETPP_TEST_LIBS})
# we can just link the legacy here because we usually don't need to go through a debugger on these tests
if (NOT GLMNETPP_MOCK_LEGACY)
    target_link_libraries(elnet_driver_unittest PRIVATE legacy)
endif()
add_test(elnet_driver_unittest elnet_driver_unittest)
