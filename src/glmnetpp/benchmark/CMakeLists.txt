# All macro tests
set(
    BENCHMARKS
    #lasso_stress_benchmark
    gaussian_cov_benchmark
    gaussian_naive_benchmark
    gaussian_benchmark
)

foreach( benchmark ${BENCHMARKS} )
    add_executable(${benchmark} ${CMAKE_CURRENT_SOURCE_DIR}/${benchmark}.cpp)
    target_compile_options(${benchmark} PRIVATE -std=c++17)
    target_include_directories(${benchmark} PRIVATE
        ${PROJECT_SOURCE_DIR}/test
        ${GLMNETPP_SOURCEDIR})
    target_link_libraries(${benchmark}
        benchmark::benchmark_main
        ${PROJECT_NAME}
        legacy
        Eigen3::Eigen)
    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
		target_link_libraries(${benchmark} pthread)
    endif()
endforeach()
