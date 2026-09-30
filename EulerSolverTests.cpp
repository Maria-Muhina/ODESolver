#include <gtest/gtest.h>
#include <cmath>

#include "EulerSolver.h"
#include "ExponentialGrowth.h"
#include "DoubleGrowth.h"
#include "TimeDependentSystem.h"

TEST(EulerSolverTest, SolvesExponentialGrowth) {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_NEAR(result[0].t, 0.0, 1e-10);
    EXPECT_NEAR(result[0].x, 1.0, 1e-10);

    EXPECT_NEAR(result[1].t, 0.1, 1e-10);
    EXPECT_NEAR(result[1].x, 1.1, 1e-10);

    EXPECT_NEAR(result[2].t, 0.2, 1e-10);
    EXPECT_NEAR(result[2].x, 1.21, 1e-10);

    EXPECT_NEAR(result[3].t, 0.3, 1e-10);
    EXPECT_NEAR(result[3].x, 1.331, 1e-10);
}

TEST(EulerSolverTest, StopsExactlyAtEndTime) {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 1.0, 0.3);

    EXPECT_EQ(result.size(), 5);

    EXPECT_NEAR(result[0].t, 0.0, 1e-10);
    EXPECT_NEAR(result[0].x, 1.0, 1e-10);

    EXPECT_NEAR(result[1].t, 0.3, 1e-10);
    EXPECT_NEAR(result[1].x, 1.3, 1e-10);

    EXPECT_NEAR(result[2].t, 0.6, 1e-10);
    EXPECT_NEAR(result[2].x, 1.69, 1e-10);

    EXPECT_NEAR(result[3].t, 0.9, 1e-10);
    EXPECT_NEAR(result[3].x, 2.197, 1e-10);

    EXPECT_NEAR(result[4].t, 1.0, 1e-10);
    EXPECT_NEAR(result[4].x, 2.4167, 1e-10);
}

TEST(EulerSolverTest, ErrorDecreasesWithSmallerStep) {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;

    double exact = std::exp(1.0);


    auto result_02 = solver.solve(system, 0.0, 1.0, 1.0, 0.2);
    double numerical_02 = result_02.back().x;
    double error_02 = std::abs(numerical_02 - exact);

    auto result_01 = solver.solve(system, 0.0, 1.0, 1.0, 0.1);
    double numerical_01 = result_01.back().x;
    double error_01 = std::abs(numerical_01 - exact);

    auto result_005 = solver.solve(system, 0.0, 1.0, 1.0, 0.05);
    double numerical_005 = result_005.back().x;
    double error_005 = std::abs(numerical_005 - exact);

    EXPECT_LT(error_01, error_02);
    EXPECT_LT(error_005, error_01);

    double p_1 = std::log2(error_02 / error_01);
    double p_2 = std::log2(error_01 / error_005);

    // Smaller h gives a more accurate estimate of the asymptotic order.
    EXPECT_NEAR(p_1, 1.0, 0.2);
    EXPECT_NEAR(p_2, 1.0, 0.1);
}

TEST(EulerSolverTest, SolvesExponentialGrowth_changed_x0) {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 2.0, 0.5, 0.1);

    EXPECT_EQ(result.size(), 6);

    EXPECT_NEAR(result[0].t, 0.0, 1e-10);
    EXPECT_NEAR(result[0].x, 2.0, 1e-10);

    EXPECT_NEAR(result[1].t, 0.1, 1e-10);
    EXPECT_NEAR(result[1].x, 2.2, 1e-10);

    EXPECT_NEAR(result[2].t, 0.2, 1e-10);
    EXPECT_NEAR(result[2].x, 2.42, 1e-10);

    EXPECT_NEAR(result[3].t, 0.3, 1e-10);
    EXPECT_NEAR(result[3].x, 2.662, 1e-10);

    EXPECT_NEAR(result[4].t, 0.4, 1e-10);
    EXPECT_NEAR(result[4].x, 2.9282, 1e-10);

    EXPECT_NEAR(result[5].t, 0.5, 1e-10);
    EXPECT_NEAR(result[5].x, 3.22102, 1e-10);
}

TEST(EulerSolverTest, SolvesDoubleGrowth) {
    ode::EulerSolver solver;
    ode::DoubleGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 1.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_NEAR(result[0].t, 0.0, 1e-10);
    EXPECT_NEAR(result[0].x, 1.0, 1e-10);

    EXPECT_NEAR(result[1].t, 0.1, 1e-10);
    EXPECT_NEAR(result[1].x, 1.2, 1e-10);

    EXPECT_NEAR(result[2].t, 0.2, 1e-10);
    EXPECT_NEAR(result[2].x, 1.44, 1e-10);

    EXPECT_NEAR(result[3].t, 0.3, 1e-10);
    EXPECT_NEAR(result[3].x, 1.728, 1e-10);
}

TEST(EulerSolverTest, SolvesTimeDependentSystem) {
    ode::EulerSolver solver;
    ode::TimeDependentSystem system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, 0.0, 0.3, 0.1);

    EXPECT_EQ(result.size(), 4);

    EXPECT_NEAR(result[0].t, 0.0, 1e-10);
    EXPECT_NEAR(result[0].x, 0.0, 1e-10);

    EXPECT_NEAR(result[1].t, 0.1, 1e-10);
    EXPECT_NEAR(result[1].x, 0.0, 1e-10);

    EXPECT_NEAR(result[2].t, 0.2, 1e-10);
    EXPECT_NEAR(result[2].x, 0.01, 1e-10);

    EXPECT_NEAR(result[3].t, 0.3, 1e-10);
    EXPECT_NEAR(result[3].x, 0.03, 1e-10);
}

TEST(EulerSolverTest, SolvesNegative_x0) {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result =
    solver.solve(system, 0.0, -1.0, 0.2, 0.1);

    EXPECT_EQ(result.size(), 3);

    EXPECT_NEAR(result[0].t, 0.0, 1e-10);
    EXPECT_NEAR(result[0].x, -1.0, 1e-10);

    EXPECT_NEAR(result[1].t, 0.1, 1e-10);
    EXPECT_NEAR(result[1].x, -1.1, 1e-10);

    EXPECT_NEAR(result[2].t, 0.2, 1e-10);
    EXPECT_NEAR(result[2].x, -1.21, 1e-10);
}