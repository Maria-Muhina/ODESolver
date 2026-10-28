#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include "OdeTypes.h"
#include "AdaptiveRK4Solver.h"
#include "ExponentialGrowth.h"
#include "HarmonicOscillator.h"

TEST(AdaptiveRK4SolverTest, ExponentialGrowth) {
    ode::AdaptiveRK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result =
        solver.solve(system, 0.0, x0, 1.0, 0.2, 1e-6);

    ASSERT_GE(result.size(), 2u);
    EXPECT_NEAR(result[0].x(0), 1.0, 1e-12);
    EXPECT_NEAR(result[0].t, 0.0, 1e-14);
    EXPECT_NEAR(result.back().x(0), std::exp(1.0), 1e-5);
    EXPECT_NEAR(result.back().t, 1, 1e-14);

    for (std::size_t i = 1; i < result.size(); i++) {
        EXPECT_GT(result[i].t, result[i-1].t);
    }
}

TEST(AdaptiveRK4SolverTest, RejectsInvalidParameters) {
    ode::AdaptiveRK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    EXPECT_THROW(
        solver.solve(system, 0.0, x0, 1.0, 0.2, 0.0),
        std::invalid_argument
    );

    EXPECT_THROW(
        solver.solve(system, 0.0, x0, 1.0, 0.2, -1e-6),
        std::invalid_argument
    );

    EXPECT_THROW(
        solver.solve(system, 0.0, x0, 1.0, 0.2, std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument
    );

    EXPECT_THROW(
        solver.solve(system, 0.0, x0, 1.0, 0.2, std::numeric_limits<double>::infinity()),
        std::invalid_argument
    );

    ode::Trajectory result = solver.solve(system, 1.0, x0, 1.0, 0.2, 0.1);

    ASSERT_EQ(result.size(), 1);
    EXPECT_DOUBLE_EQ(result[0].t, 1.0);
    EXPECT_DOUBLE_EQ(result[0].x(0), 1.0);
}

TEST(AdaptiveRK4SolverTest, RejectsLargeInitialStep) {
    ode::AdaptiveRK4Solver solver;
    ode::ExponentialGrowth system;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result =
        solver.solve(system, 0.0, x0, 1.0, 1.0, 1e-9);

    ASSERT_GT(result.size(), 2u);
    EXPECT_GT(1.0, result[1].t);
    EXPECT_NEAR(result.back().t, 1.0, 1e-14);

    for (std::size_t i = 1; i < result.size(); i++) {
        EXPECT_GT(result[i].t, result[i-1].t);
    }
}

class ZeroSystem : public ode::OdeSystem {
public:
    void rhs(
        double /*t*/,
        const ode::State& x,
        ode::State& dxdt
    ) const override {
        dxdt = ode::State::Zero(x.size());
    }
};

TEST(AdaptiveRK4SolverTest, HandlesZeroError) {
    ode::AdaptiveRK4Solver solver;
    ZeroSystem system;

    ode::State x0(2);
    x0 << 1.0, 2.0;

    ode::Trajectory result =
        solver.solve(system, 0.0, x0, 1.0, 0.1, 1e-6);

    ASSERT_GE(result.size(), 3u);
    EXPECT_NEAR(result[1].t, 0.1, 1e-14);
    EXPECT_NEAR(result[2].t, 0.6, 1e-14);
    EXPECT_NEAR(result.back().t, 1.0, 1e-14);
    EXPECT_NEAR(result.back().x(0), 1.0, 1e-14);
    EXPECT_NEAR(result.back().x(1), 2.0, 1e-14);

    for (std::size_t i = 1; i < result.size(); i++) {
        EXPECT_GT(result[i].t, result[i-1].t);
    }
}

TEST(AdaptiveRK4SolverTest, SolvesHarmonicOscillator) {
    ode::AdaptiveRK4Solver solver;
    ode::HarmonicOscillator system;

    ode::State x0(2);
    x0 << 1.0, 0.0;

    ode::Trajectory result =
        solver.solve(system, 0.0, x0, 10.0, 0.5, 1e-7);

    ASSERT_GE(result.size(), 2u);
    ASSERT_EQ(result.back().x.size(), 2);
    EXPECT_NEAR(result.back().t, 10.0, 1e-14);
    EXPECT_NEAR(result.back().x(0), std::cos(10.0), 1e-4);
    EXPECT_NEAR(result.back().x(1), -std::sin(10.0), 1e-4);

    for (std::size_t i = 1; i < result.size(); i++) {
        EXPECT_GT(result[i].t, result[i-1].t);
    }
}

TEST(AdaptiveRK4SolverTest, RejectsNonAdvancingTime) {
    ode::AdaptiveRK4Solver solver;
    ZeroSystem system;

    ode::State x0(2);
    x0 << 1.0, 2.0;

    double t0 = 1e16;
    double t_end = t0 + 4.0;
    double h_initial = 1.0;

    EXPECT_THROW(
        solver.solve(system, t0, x0, t_end, h_initial, 1e-6),
        std::runtime_error
    );
}
