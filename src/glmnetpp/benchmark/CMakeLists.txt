# All macro tests
set(
    BENCHMARKS
    poisson_benchmark
    binomial_benchmark
    binomial_two_class_benchmark
    gaussian_cov_benchmark
    gaussian_naive_benchmark
    gaussian_benchmark
)

foreach( benchmark ${BENCHMARKS} )
	set(BENCHMARK_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${benchmark}.cpp)
    if (NOT GLMNETPP_MOCK_LEGACY) 
        set(BENCHMARK_SOURCES 
        ${BENCHMARK_SOURCES}
        ${GLMNETPP_SOURCEDIR}/legacy/glmnet5dpclean.f
        ${GLMNETPP_SOURCEDIR}/legacy/pb.c   # only needed for fortran code to link to this
            ) 
    endif()

    add_executable(${benchmark} ${BENCHMARK_SOURCES})
    set_property(TARGET ${benchmark} PROPERTY CXX_STANDARD 17)
    target_include_directories(${benchmark} PRIVATE
        ${PROJECT_SOURCE_DIR}/test
        ${GLMNETPP_SOURCEDIR})
    target_link_libraries(${benchmark}
        benchmark::benchmark_main
        ${PROJECT_NAME}
        Eigen3::Eigen)
	# This is if we can't compile Fortran code,
	# but we still want to run our C++ method.
	# We enable the mocked version of the legacy Fortran code.
	if (GLMNETPP_MOCK_LEGACY)
		target_compile_definitions(${benchmark} PRIVATE
		-DGLMNETPP_MOCK_LEGACY)
    else ()
        target_include_directories(${benchmark} PRIVATE ${R_INCLUDE_DIR})
        target_link_libraries(${benchmark} ${RLIB})
	endif()
    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
		target_link_libraries(${benchmark} pthread)
    endif()
endforeach()
