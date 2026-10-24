#include <cassert>
#include <cmath>
#include "OdeSystem.h"
#include "ExponentialGrowth.h"
#include "DoubleGrowth.h"
#include "EulerSolver.h"

void test_exponential_growth() {
    ode::ExponentialGrowth system;

    ode::State x(1);
    x << 1.0;
    ode::State expected(1);
    expected << 1.0;

    ode::check_system(system, 0.0, x, expected);
    x(0) = 3.0;
    expected(0) = 3.0;
    ode::check_system(system, 2.0, x, expected);
    x(0) = -5.0;
    expected(0) = -5.0;
    ode::check_system(system, 10.0, x, expected);
}

void test_double_growth() {
    ode::DoubleGrowth system;

    ode::State x(1);
    x << 1.0;
    ode::State expected(1);
    expected << 2.0;

    ode::check_system(system, 0.0, x, expected);
    x(0) = 3.0;
    expected(0) = 6.0;
    ode::check_system(system, 2.0, x, expected);
    x(0) = -5.0;
    expected(0) = -10.0;
    ode::check_system(system, 10.0, x, expected);
}

void test_check_system() {
    ode::ExponentialGrowth exponential;
    ode::DoubleGrowth double_growth;

    ode::State x(1);
    x << 1.0;
    ode::State expected(1);
    expected << 1.0;

    ode::check_system(exponential, 0.0, x, expected);
    expected(0) = 2.0;
    ode::check_system(double_growth, 0.0, x, expected);
}

void test_euler_solver() {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;
    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory result = 
        solver.solve(system, 0.0, x0, 0.3, 0.1);

    assert(std::abs(result[0].t - 0.0) < 1e-10);
    assert(std::abs(result[0].x(0) - 1.0) < 1e-10);
    
    assert(std::abs(result[1].t - 0.1) < 1e-10);
    assert(std::abs(result[1].x(0) - 1.1) < 1e-10);

    assert(std::abs(result[2].t - 0.2) < 1e-10);
    assert(std::abs(result[2].x(0) - 1.21) < 1e-10);

    assert(std::abs(result[3].t - 0.3) < 1e-10);
    assert(std::abs(result[3].x(0) - 1.331) < 1e-10);

    assert(result.size() == 4);

}

int main() {

    test_exponential_growth();
    test_double_growth();
    test_check_system();
    test_euler_solver();

    return 0;
}