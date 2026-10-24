#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "OdeTypes.h"
#include "RK4Solver.h"
#include "ExponentialGrowth.h"
#include "DoubleGrowth.h"
#include "TimeDependentSystem.h"
#include "HarmonicOscillator.h"

TEST(RK4SolverTest, MethodTest) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result =
    solver.solve(system, 0.0, x0, 0.1, 0.1);

    EXPECT_EQ(result.size(), 2);
    EXPECT_NEAR(result[1].x(0), 1.1051708333333333, 1e-12);
}

TEST(RK4SolverTest, ExactSolutionCompareTest) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result =
    solver.solve(system, 0.0, x0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_DOUBLE_EQ(result[0].x(0), 1.0);
    EXPECT_DOUBLE_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x(0), std::exp(0.1), 1e-6);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x(0), std::exp(0.2), 1e-6);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x(0), std::exp(0.3), 1e-6);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, ExactSolutionCompareTest_changed_x0) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 2.0;

    ode::Trajectory result =
    solver.solve(system, 0.0, x0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_DOUBLE_EQ(result[0].x(0), 2.0);
    EXPECT_DOUBLE_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x(0), 2.0 * std::exp(0.1), 1e-6);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x(0), 2.0 * std::exp(0.2), 1e-6);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x(0), 2.0 * std::exp(0.3), 1e-6);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, SolvesDoubleGrowth) {
    ode::RK4Solver solver;
    ode::DoubleGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result =
    solver.solve(system, 0.0, x0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_DOUBLE_EQ(result[0].x(0), 1.0);
    EXPECT_DOUBLE_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x(0), std::exp(2 * 0.1), 1e-4);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x(0), std::exp(2 * 0.2), 1e-4);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x(0), std::exp(2 * 0.3), 1e-4);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, SolvesTimeDependentSystem) {
    ode::RK4Solver solver;
    ode::TimeDependentSystem system;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result =
    solver.solve(system, 0.0, x0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_DOUBLE_EQ(result[0].x(0), 1.0);
    EXPECT_DOUBLE_EQ(result[0].t, 0.0);

    EXPECT_NEAR(result[1].x(0), 1.005, 1e-14);
    EXPECT_NEAR(result[1].t, 0.1, 1e-7);

    EXPECT_NEAR(result[2].x(0), 1.02, 1e-14);
    EXPECT_NEAR(result[2].t, 0.2, 1e-7);

    EXPECT_NEAR(result[3].x(0), 1.045, 1e-14);
    EXPECT_NEAR(result[3].t, 0.3, 1e-7);
}

TEST(RK4SolverTest, ErrorDecreasesWithSmallerStep) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    double exact = std::exp(0.3);

    ode::Trajectory result_02 =
        solver.solve(system, 0.0, x0, 0.3, 0.2);
    double numerical_02 = result_02.back().x(0);
    double error_02 = std::abs(exact - numerical_02);
    ASSERT_LT(error_02, 1e-5);

    ode::Trajectory result_01 =
        solver.solve(system, 0.0, x0, 0.3, 0.1);
    double numerical_01 = result_01.back().x(0);
    double error_01 = std::abs(exact - numerical_01);
    ASSERT_LT(error_01, 1e-6);

    ode::Trajectory result_005 =
        solver.solve(system, 0.0, x0, 0.3, 0.05);
    double numerical_05 = result_005.back().x(0);
    double error_005 = std::abs(exact - numerical_05);
    ASSERT_LT(error_005, 1e-6);

    ode::Trajectory result_0025 =
        solver.solve(system, 0.0, x0, 0.3, 0.025);
    double numerical_025 = result_0025.back().x(0);
    double error_0025 = std::abs(exact - numerical_025);
    ASSERT_LT(error_0025, 1e-6);

    double p1 = std::log(error_02 / error_01)/std::log(2.0);
    EXPECT_NEAR(p1, 4.0, 0.7);

    double p2 = std::log(error_01 / error_005)/std::log(2.0);
    EXPECT_NEAR(p2, 4.0, 0.1);

    double p3 = std::log(error_005 / error_0025)/std::log(2.0);
    EXPECT_NEAR(p3, 4.0, 0.05);
}

TEST(RK4SolverTest, SolvesHarmonicOscillator) {
    ode::RK4Solver solver;
    ode::HarmonicOscillator system;

    ode::State x0(2);
    x0 << 1.0, 0.0;

    ode::Trajectory result =
    solver.solve(system, 0.0, x0, 1.0, 0.1);

    EXPECT_EQ(result.size(), 11);

    EXPECT_EQ(result.back().x.size(), 2);
    EXPECT_NEAR(result.back().x(0), std::cos(1.0), 1e-4);
    EXPECT_NEAR(result.back().x(1), -std::sin(1.0), 1e-4);
    EXPECT_NEAR(result.back().t, 1.0, 1e-7);
}

TEST(RK4SolverTest, RejectsInvalidTimeParameters) {
    ode::RK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    EXPECT_THROW(
        solver.solve(system, 0.0, x0, 1.0, 0.0),
        std::invalid_argument
    );

    EXPECT_THROW(
        solver.solve(system, 0.0, x0, 1.0, -0.1),
        std::invalid_argument
    );

    EXPECT_THROW(
        solver.solve(system, 2.0, x0, 1.0, 0.1),
        std::invalid_argument
    );

    EXPECT_THROW(
        solver.solve(
            system, 0.0, x0, 1.0,
            std::numeric_limits<double>::quiet_NaN()
        ),
        std::invalid_argument
    );

    ode::Trajectory result = solver.solve(system, 1.0, x0, 1.0, 0.1);

    ASSERT_EQ(result.size(), 1);
    EXPECT_DOUBLE_EQ(result[0].t, 1.0);
    EXPECT_DOUBLE_EQ(result[0].x(0), 1.0);
}
