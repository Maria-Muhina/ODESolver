#include <gtest/gtest.h>
#include <cmath>

#include "RK4Solver.h"
#include "ExponentialGrowth.h"

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