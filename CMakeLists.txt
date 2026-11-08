cmake_minimum_required(VERSION 3.20)
project(ASC_bla LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include_directories(src)
include_directories(HPC/src)

enable_testing()

# Release builds define NDEBUG, which removes the element bounds checks.
add_executable(test_bla tests/test_bla.cc)
add_test(NAME test_bla COMMAND test_bla)

# The same tests with bounds checks enabled, whatever the build type.
add_executable(test_bla_checked tests/test_bla.cc)
target_compile_options(test_bla_checked PRIVATE -UNDEBUG)
add_test(NAME test_bla_checked COMMAND test_bla_checked)

# Demonstration programs that print their results.
add_executable(test_vector tests/test_vector.cc)
add_executable(test_matrix tests/test_matrix.cc)

set(BLA_SIZEOF_INTEGER 4)
find_package(LAPACK QUIET)
if(LAPACK_FOUND)
  add_executable(test_lapack tests/test_lapack.cc)
  target_link_libraries(test_lapack PUBLIC ${LAPACK_LIBRARIES})
  add_executable(test_lapack_checks tests/test_lapack_checks.cc)
  target_link_libraries(test_lapack_checks PUBLIC ${LAPACK_LIBRARIES})
  add_test(NAME test_lapack_checks COMMAND test_lapack_checks)
else()
  message(STATUS "LAPACK not found: skipping the LAPACK tests")
endif()

# Python module, built when pybind11 is available.
find_package(Python COMPONENTS Interpreter Development.Module QUIET)
if(Python_FOUND)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR ERROR_QUIET)
  list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
endif()
find_package(pybind11 CONFIG QUIET)
if(pybind11_FOUND)
  pybind11_add_module(bla src/bind_bla.cpp)
  install(TARGETS bla DESTINATION ASCsoft)
  install(FILES src/vector.h DESTINATION ASCsoft/include)
else()
  message(STATUS "pybind11 not found: skipping the Python module")
endif()
