#include <gtest/gtest.h>

#include "EulerSolver.h"
#include "ExponentialGrowth.h"

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