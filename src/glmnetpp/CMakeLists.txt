cmake_minimum_required(VERSION 3.7)

project("glmnetpp" VERSION 4.2.0
                 DESCRIPTION "A C++ implementation of elastic net solver.")

option(GLMNETPP_ENABLE_TEST         "Enable unit tests to be built."    ON)
option(GLMNETPP_ENABLE_BENCHMARK    "Enable benchmarks to be built."    OFF)
option(GLMNETPP_ENABLE_COVERAGE     "Build glmnetpp with coverage"      OFF)
option(GLMNETPP_MOCK_LEGACY        "Build glmnetpp with mock versions of the fortran."      OFF)

# Define a set of strict build flags
set(GLMNETPP_STRICT_WARNINGS -g -Wall)
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(GLMNETPP_STRICT_WARNINGS
        ${GLMNETPP_STRICT_WARNINGS}
        -Wextra -Werror -pedantic)
endif()

# Stoopid hack on windows
if (WIN32)
	SET(CMAKE_FIND_LIBRARY_PREFIXES "")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
endif()

# Dependency on Eigen
find_package(Eigen3 3.3 CONFIG REQUIRED NO_MODULE)
message(STATUS "Eigen3 found at ${EIGEN3_INCLUDE_DIR}")

# Set include dirs
set(GLMNETPP_INCLUDEDIR "${${PROJECT_NAME}_SOURCE_DIR}/include")
set(GLMNETPP_SOURCEDIR "${${PROJECT_NAME}_SOURCE_DIR}/src")

# Add this library as shared library.
add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} SYSTEM INTERFACE
            $<BUILD_INTERFACE:${GLMNETPP_INCLUDEDIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Set C++14 standard for project target
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_14)

# Set install destinations
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Create GlmnetppConfigVersion.cmake which contains current project version
# This is supposed to help with (major) version compatibility.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

# Build source library
if (GLMNETPP_ENABLE_TEST OR GLMNETPP_ENABLE_BENCHMARK)
    if (NOT GLMNETPP_MOCK_LEGACY)
		# MUST be called at the top-level CMakeLists.txt
		enable_language(Fortran)
    endif()
endif()
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# Automate the choosing of config
# if CMAKE_BUILD_TYPE not defined
if (NOT CMAKE_BUILD_TYPE)
    # if binary directory ends with "release", use release mode
    if (${PROJECT_BINARY_DIR} MATCHES "release$")
        set(CMAKE_BUILD_TYPE RELEASE)
    # otherwise, use debug mode
    else()
        set(CMAKE_BUILD_TYPE DEBUG)
    endif()
endif()
message(STATUS "Compiling in ${CMAKE_BUILD_TYPE} mode")

if (GLMNETPP_ENABLE_TEST)
    # This will perform memcheck
    include(CTest)
    enable_testing()

    # Find googletest
    find_package(GTest 1.11 CONFIG REQUIRED)

    # add test subdirectory
    add_subdirectory(${PROJECT_SOURCE_DIR}/test ${PROJECT_BINARY_DIR}/test)
endif()

if (GLMNETPP_ENABLE_BENCHMARK)
    find_package(benchmark REQUIRED
        HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../../benchmark/build)

    # add benchmark subdirectory
    add_subdirectory(${PROJECT_SOURCE_DIR}/benchmark ${PROJECT_BINARY_DIR}/benchmark)
endif()
