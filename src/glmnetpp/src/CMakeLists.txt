########################################################################
# Legacy Library for compatability tests/benchmarks
########################################################################
if (NOT GLMNETPP_MOCK_LEGACY)
    add_library(legacy
        ${CMAKE_CURRENT_SOURCE_DIR}/legacy/wls.f
        ${CMAKE_CURRENT_SOURCE_DIR}/legacy/glmnet5dpclean.f
        ${CMAKE_CURRENT_SOURCE_DIR}/legacy/pb.c
        )
    if (NOT R_INCLUDE_DIR)
        message(FATAL_ERROR
            "R include directory must be provided by the CMake variable R_INCLUDE_DIR. "
            "Run Sys.getenv('R_INCLUDE_DIR') in R interpreter and provide this path. ")
    endif()
    if (NOT R_LIB_DIR)
        message(FATAL_ERROR
            "R lib directory must be provided by the CMake variable R_LIB_DIR. "
            "Run Sys.getenv('R_HOME') in R interpreter, append /lib, "
            "and provide the resulting path. ")
    endif()
    target_compile_options(legacy PRIVATE -Wall -Wextra)
    target_include_directories(legacy PRIVATE ${R_INCLUDE_DIR})
    find_library(RLIB R REQUIRED HINTS ${R_LIB_DIR})
    target_link_libraries(legacy PRIVATE ${RLIB})
endif()
