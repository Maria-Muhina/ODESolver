cmake_minimum_required(VERSION 3.16)

project(StateProject)

add_executable(main
    main.cpp
    OdeSystem.cpp
    ExponentialGrowth.cpp
    DoubleGrowth.cpp
    Counter.cpp
    EulerSolver.cpp
)