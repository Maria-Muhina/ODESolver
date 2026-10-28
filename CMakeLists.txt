cmake_minimum_required(VERSION 3.16)

project(ODESolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

enable_testing()

include(FetchContent)

FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
)

FetchContent_MakeAvailable(eigen)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

FetchContent_MakeAvailable(googletest)

add_executable(main
    main.cpp
    OdeSystem.cpp
    ExponentialGrowth.cpp
    DoubleGrowth.cpp
    EulerSolver.cpp
    TimeDependentSystem.cpp
    RK4Solver.cpp
    OdeValidation.cpp
    RK4Step.cpp
    AdaptiveRK4Solver.cpp
)

add_executable(ODESolverTests
    EulerSolverTests.cpp
    RK4SolverTests.cpp
    AdaptiveRK4SolverTests.cpp
    OdeSystem.cpp
    ExponentialGrowth.cpp
    DoubleGrowth.cpp
    EulerSolver.cpp
    TimeDependentSystem.cpp
    RK4Solver.cpp
    HarmonicOscillator.cpp
    OdeValidation.cpp
    RK4Step.cpp
    AdaptiveRK4Solver.cpp
)

add_executable(Eigen_test Eigen_test.cpp)

target_link_libraries(Eigen_test PRIVATE Eigen3::Eigen)

target_link_libraries(ODESolverTests
    PRIVATE
    GTest::gtest_main
    Eigen3::Eigen
)

target_link_libraries(main PRIVATE Eigen3::Eigen)

include(GoogleTest)

gtest_discover_tests(ODESolverTests)