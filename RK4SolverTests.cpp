#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "RK4Solver.h"
#include "ExponentialGrowth.h"
#include "DoubleGrowth.h"
#include "TimeDependentSystem.h"

TEST(RK4SolverTest, MethodTest) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 0.1, 0.1);

    EXPECT_EQ(result.size(), 2);
    EXPECT_NEAR(result[1].x, 1.1051708333333333, 1e-12);
}

TEST(RK4SolverTest, ExactSolutionCompareTest) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_EQ(result[0].x, 1.0);
    EXPECT_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x, std::exp(0.1), 1e-6);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x, std::exp(0.2), 1e-6);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x, std::exp(0.3), 1e-6);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, ExactSolutionCompareTest_changed_x0) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 2.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_EQ(result[0].x, 2.0);
    EXPECT_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x, 2.0 * std::exp(0.1), 1e-6);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x, 2.0 * std::exp(0.2), 1e-6);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x, 2.0 * std::exp(0.3), 1e-6);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, SolvesDoubleGrowth) {
    ode::RK4Solver solver;
    ode::DoubleGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_EQ(result[0].x, 1.0);
    EXPECT_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x, std::exp(2 * 0.1), 1e-4);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x, std::exp(2 * 0.2), 1e-4);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x, std::exp(2 * 0.3), 1e-4);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, SolvesTimeDependentSystem) {
    ode::RK4Solver solver;
    ode::TimeDependentSystem system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_EQ(result[0].x, 1.0);
    EXPECT_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x, 1.005, 1e-14);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x, 1.02, 1e-14);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x, 1.045, 1e-14);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, ErrorDecreasesWithSmallerStep) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    double exact = std::exp(0.3);

    std::vector<ode::SolutionPoint> result_02 =
        solver.solve(system, 0.0, 1.0, 0.3, 0.2);
    double numerical_02 = result_02.back().x;
    double error_02 = std::abs(exact - numerical_02);
    ASSERT_LT(error_02, 1e-5);

    std::vector<ode::SolutionPoint> result_01 =
        solver.solve(system, 0.0, 1.0, 0.3, 0.1);
    double numerical_01 = result_01.back().x;
    double error_01 = std::abs(exact - numerical_01);
    ASSERT_LT(error_01, 1e-6);

    std::vector<ode::SolutionPoint> result_005 =
        solver.solve(system, 0.0, 1.0, 0.3, 0.05);
    double numerical_05 = result_005.back().x;
    double error_005 = std::abs(exact - numerical_05);
    ASSERT_LT(error_005, 1e-6);

    std::vector<ode::SolutionPoint> result_0025 =
        solver.solve(system, 0.0, 1.0, 0.3, 0.025);
    double numerical_025 = result_0025.back().x;
    double error_0025 = std::abs(exact - numerical_025);
    ASSERT_LT(error_0025, 1e-6);

    double p1 = std::log(error_02 / error_01)/std::log(2.0);
    EXPECT_NEAR(p1, 4.0, 0.7);

    double p2 = std::log(error_01 / error_005)/std::log(2.0);
    EXPECT_NEAR(p2, 4.0, 0.1);

    double p3 = std::log(error_005 / error_0025)/std::log(2.0);
    EXPECT_NEAR(p3, 4.0, 0.05);

    EXPECT_LT(p1, p2);
    EXPECT_LT(p2, p3);
}
