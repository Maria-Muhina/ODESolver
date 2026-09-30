cmake_minimum_required(VERSION 3.16)

project(ODESolver)

enable_testing()

include(FetchContent)

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
)

add_executable(EulerSolverTests
    EulerSolverTests.cpp
    OdeSystem.cpp
    ExponentialGrowth.cpp
    DoubleGrowth.cpp
    EulerSolver.cpp
    TimeDependentSystem.cpp
)

target_link_libraries(EulerSolverTests
    PRIVATE
    GTest::gtest_main
)

include(GoogleTest)

gtest_discover_tests(EulerSolverTests)