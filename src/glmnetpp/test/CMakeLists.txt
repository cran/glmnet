# Try to find package Threads for the use of pthreads
unset(GLMNETPP_HAS_PTHREAD)
find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
    set(GLMNETPP_HAS_PTHREAD ON)
endif()

if (DEFINED GLMNETPP_HAS_PTHREAD)
    set(GLMNETPP_HAS_PTHREAD_MACRO "-DGLMNETPP_HAS_PTHREAD")
endif()

########################################################################
# Utility TEST
########################################################################

add_executable(utility_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/util/type_traits_unittest.cpp
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(utility_unittest PRIVATE -Werror -Wextra -Wpedantic)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(utility_unittest PRIVATE -fopenmp)
    target_link_libraries(utility_unittest -fopenmp)
endif()
target_compile_options(utility_unittest PRIVATE -g -Wall)
target_include_directories(utility_unittest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
if (GLMNETPP_ENABLE_COVERAGE)
    target_link_libraries(utility_unittest gcov)
endif()
target_link_libraries(utility_unittest
    ${PROJECT_NAME}
    GTest::gtest_main
    Eigen3::Eigen)
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(utility_unittest pthread)
endif()
add_test(utility_unittest utility_unittest)

########################################################################
# Translation TEST
########################################################################

add_executable(translation_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet1_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/elnet2_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/spelnet_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/spelnet1_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/translation/spelnet2_unittest.cpp
    # Explicitly add Fortran source for debugging purposes
    ${GLMNETPP_SOURCEDIR}/legacy/glmnet5dpclean.f
    ${GLMNETPP_SOURCEDIR}/legacy/pb.c   # only needed for fortran code to link to this
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # DO NOT put -Werror because Fortran gives a bunch of warnings that we unfortunately cannot change.
    target_compile_options(translation_unittest PRIVATE -Wextra -Wpedantic)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(translation_unittest PRIVATE -fopenmp)
    target_link_libraries(translation_unittest -fopenmp)
endif()
target_compile_options(translation_unittest PRIVATE -g -Wall)
target_compile_definitions(translation_unittest PRIVATE
    ${GLMNETPP_HAS_PTHREAD_MACRO})
# This is if we can't compile Fortran code,
# but we still want to run our C++ method.
# We enable the mocked version of the legacy Fortran code.
if (GLMNETPP_MOCK_LEGACY)
    target_compile_definitions(translation_unittest PRIVATE
    -DGLMNETPP_MOCK_LEGACY)
endif()
target_include_directories(translation_unittest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${R_INCLUDE_DIR}
    ${GLMNETPP_SOURCEDIR})
if (GLMNETPP_ENABLE_COVERAGE)
    target_link_libraries(translation_unittest gcov)
endif()
target_link_libraries(translation_unittest
    ${PROJECT_NAME}
    GTest::gtest_main
    Eigen3::Eigen)
if (NOT GLMNETPP_MOCK_LEGACY)
    target_link_libraries(translation_unittest ${RLIB})
endif()
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(translation_unittest pthread)
endif()
add_test(translation_unittest translation_unittest)

########################################################################
# Elnet Path TEST
########################################################################

add_executable(elnet_path_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/gaussian_cov_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/gaussian_naive_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_gaussian_cov_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_path/sp_gaussian_naive_unittest.cpp
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(elnet_path_unittest PRIVATE -Werror -Wextra -Wpedantic)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(elnet_path_unittest PRIVATE -fopenmp)
    target_link_libraries(elnet_path_unittest -fopenmp)
endif()
target_compile_options(elnet_path_unittest PRIVATE -g -Wall)
target_compile_definitions(elnet_path_unittest PRIVATE
    ${GLMNETPP_HAS_PTHREAD_MACRO})
target_include_directories(elnet_path_unittest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GLMNETPP_SOURCEDIR})
if (GLMNETPP_ENABLE_COVERAGE)
    target_link_libraries(elnet_path_unittest gcov)
endif()
target_link_libraries(elnet_path_unittest
    ${PROJECT_NAME}
    GTest::gtest_main
    Eigen3::Eigen)
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(elnet_path_unittest pthread)
endif()
add_test(elnet_path_unittest elnet_path_unittest)

########################################################################
# Elnet Driver TEST (Integration Test)
########################################################################

add_executable(elnet_driver_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/elnet_driver/gaussian_unittest.cpp
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(elnet_driver_unittest PRIVATE -Werror -Wextra -Wpedantic)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(elnet_driver_unittest PRIVATE -fopenmp)
    target_link_libraries(elnet_driver_unittest -fopenmp)
endif()
target_compile_options(elnet_driver_unittest PRIVATE -g -Wall)
target_compile_definitions(elnet_driver_unittest PRIVATE
    ${GLMNETPP_HAS_PTHREAD_MACRO})
target_include_directories(elnet_driver_unittest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GLMNETPP_SOURCEDIR})
if (GLMNETPP_ENABLE_COVERAGE)
    target_link_libraries(elnet_driver_unittest gcov)
endif()
target_link_libraries(elnet_driver_unittest
    ${PROJECT_NAME}
    GTest::gtest_main
    Eigen3::Eigen)
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(elnet_driver_unittest pthread)
endif()
add_test(elnet_driver_unittest elnet_driver_unittest)

########################################################################
# Core TEST
########################################################################

add_executable(core_unittest
    ${CMAKE_CURRENT_SOURCE_DIR}/wls_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/standardize_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chkvars_unittest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/internal_unittest.cpp
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(core_unittest PRIVATE -Werror -Wextra -Wpedantic)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(core_unittest PRIVATE -fopenmp)
    target_link_libraries(core_unittest -fopenmp)
endif()
target_compile_options(core_unittest PRIVATE -g -Wall)
target_compile_definitions(core_unittest PRIVATE ${GLMNETPP_HAS_PTHREAD_MACRO})
target_include_directories(core_unittest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GLMNETPP_SOURCEDIR})
if (GLMNETPP_ENABLE_COVERAGE)
    target_link_libraries(core_unittest gcov)
endif()
target_link_libraries(core_unittest
    ${PROJECT_NAME}
    GTest::gtest_main
    Eigen3::Eigen)
if (NOT GLMNETPP_MOCK_LEGACY)
    target_link_libraries(core_unittest legacy)
endif()
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(core_unittest pthread)
endif()
add_test(core_unittest core_unittest)
