#include <cassert>
#include <cmath>
#include "OdeSystem.h"
#include "ExponentialGrowth.h"
#include "DoubleGrowth.h"
#include "Counter.h"
#include "EulerSolver.h"

void test_exponential_growth() {
    ode::ExponentialGrowth system;

    assert(system.rhs(0.0, 1.0) == 1.0);
    assert(system.rhs(2.0, 3.0) == 3.0);
    assert(system.rhs(10.0, -5.0) == -5.0);
};

void test_double_growth() {
    ode::DoubleGrowth system;

    assert(system.rhs(0.0, 1.0) == 2.0);
    assert(system.rhs(2.0, 3.0) == 6.0);
    assert(system.rhs(10.0, -5.0) == -10.0);
};

void test_check_system() {
    ode::ExponentialGrowth exponential;
    ode::DoubleGrowth double_growth;

    ode::check_system(exponential, 0.0, 1.0, 1.0);
    ode::check_system(double_growth, 0.0, 1.0, 2.0);
}

void test_get_count() {
    Counter a;
    Counter b;
    Counter c;

    assert(Counter::get_count() == 3);
}

void test_same_count() {
    Counter a;
    Counter b;
    Counter c;

    assert(Counter::get_count() == 6);
}

void test_counter_state() {
    Counter counter;
    assert(counter.get_state() == CounterState::Inactive);

    counter.start();
    assert(counter.get_state() == CounterState::Active);

    counter.finish();
    assert(counter.get_state() == CounterState::Finished);
}

void test_euler_solver() {
    ode::EulerSolver solver;
    ode::ExponentialGrowth system;

    std::vector<ode::SolutionPoint> result = 
        solver.solve(system, 0.0, 1.0, 0.3, 0.1);

    assert(std::abs(result[0].t - 0.0) < 1e-10);
    assert(std::abs(result[0].x - 1.0) < 1e-10);
    
    assert(std::abs(result[1].t - 0.1) < 1e-10);
    assert(std::abs(result[1].x - 1.1) < 1e-10);

    assert(std::abs(result[2].t - 0.2) < 1e-10);
    assert(std::abs(result[2].x - 1.21) < 1e-10);

    assert(std::abs(result[3].t - 0.3) < 1e-10);
    assert(std::abs(result[3].x - 1.331) < 1e-10);

    assert(result.size() == 4);

}

int main() {

    test_exponential_growth();
    test_double_growth();
    test_check_system();
    test_get_count();
    test_same_count();
    test_counter_state();
    test_euler_solver();

    return 0;
}