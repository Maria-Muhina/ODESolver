#include <chrono>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "OdeSystem.h"
#include "EulerSolver.h"
#include "RK4Solver.h"
#include "ExponentialGrowth.h"
#include "AdaptiveRK4Solver.h"

class CountingSystem : public ode::OdeSystem {
public:
    explicit CountingSystem(const ode::OdeSystem& system)
        : system_(system) {}

    void rhs(
        double t,
        const ode::State& x,
        ode::State& dxdt
    ) const override {
        system_.rhs(t, x, dxdt);
        count_ ++;
    };

    std::size_t count() const {return count_;}
    void reset() const {count_ = 0;}

private:
    const ode::OdeSystem& system_;
    mutable std::size_t count_ = 0;
};

template <typename Run>
double benchmark_solver(
    const Run& run, 
    double t_end,
    int repetitions = 100
) {

    if (repetitions <= 0) {
        throw std::invalid_argument(
            "Repetitions must be positive"
        );
    }

    double checksum = 0.0;
    run();

    auto start = std::chrono::steady_clock::now();
    for (int i=0; i<repetitions; i++){
        auto result = run();
        checksum += result.back().x(0); 
    }
    auto end = std::chrono::steady_clock::now();

    std::cout << "checksum = " << checksum << "\n";
    auto result = run();
    std::cout << "steps count = " << result.size() - 1 << "\n";
    double x_end = result.back().x(0);
    double x_exact = std::exp(t_end);
    double e = std::abs(x_end - x_exact) / x_exact;
    std::cout << "error = " << e << "\n";

    std::chrono::duration<double, std::micro> elapsed = (end - start) / repetitions;
    return elapsed.count();
}

template <typename Solver>
double find_step_for_tolerance(
    const Solver& solver,
    const ode::OdeSystem& system,
    double t0,
    const ode::State& x0,
    double t_end,
    double initial_h,
    double tolerance
) {
    double h = initial_h;
    int iter = 0;

    while (h > 0 && iter < 1000) {
        auto result = solver.solve(system, t0, x0, t_end, h);
        double x_end = result.back().x(0);
        double x_exact = std::exp(t_end);
        double e = std::abs(x_end - x_exact) / x_exact;
        if (std::isfinite(e) && e <= tolerance) {
            return h;
        }
        h /= 2.0;
        iter ++;
    }
    throw std::runtime_error("Failed to find a suitable step");
}

int main() {
    ode::EulerSolver solver_euler;
    ode::RK4Solver solver_rk4;
    ode::ExponentialGrowth system;
    CountingSystem counted(system);

    ode::State x0(1);
    x0 << 1.0;
    

    double h_euler = find_step_for_tolerance(solver_euler, system, 0.0, x0, 10, 1.0, 1e-3);
    std::cout << "h_euler: " << h_euler << "\n";
    auto run_euler = [&]() {
        return solver_euler.solve(system, 0.0, x0, 10.0, h_euler);
    };
    double result_h_euler = benchmark_solver(run_euler, 10);
    std::cout
        << "Average time (eu_h_euler): " << result_h_euler << " us\n";
    std::cout << "\n";

    double h_rk4 = find_step_for_tolerance(solver_rk4, system, 0.0, x0, 10, 1.0, 1e-3);
    std::cout << "h_rk4: " << h_rk4 << "\n";
    auto run_rk4 = [&]() {
        return solver_rk4.solve(system, 0.0, x0, 10.0, h_rk4);
    };
    double result_h_rk4 = benchmark_solver(run_rk4, 10);
    std::cout << "Average time (RK4_h_rk4): " << result_h_rk4 << " us\n";

    ode::AdaptiveRK4Solver adaptive_rk4;
    auto run_adaptive_e3 = [&]() {
        return adaptive_rk4.solve(
            system, 0.0, x0, 10.0, 1.0, 1e-3
        );
    };
    auto run_adaptive_e6 = [&]() {
        return adaptive_rk4.solve(
            system, 0.0, x0, 10.0, 1.0, 1e-6
        );
    };
    auto run_adaptive_e3_counted = [&]() {
        return adaptive_rk4.solve(
            counted, 0.0, x0, 10.0, 1.0, 1e-3
        );
    };
    auto run_adaptive_e6_counted = [&]() {
        return adaptive_rk4.solve(
            counted, 0.0, x0, 10.0, 1.0, 1e-6
        );
    };

    double result_adaptive_e3 = benchmark_solver(run_adaptive_e3, 10);
    std::cout << "Average time (adaptive_e3): " << result_adaptive_e3 << " us\n";

    counted.reset();
    run_adaptive_e3_counted();
    std::cout << "rhs calls (tol = 1e-3) = " << counted.count() << "\n";

    double result_adaptive_e6 = benchmark_solver(run_adaptive_e6, 10);
    std::cout << "Average time (adaptive_e6): " << result_adaptive_e6 << " us\n";

    counted.reset();
    run_adaptive_e6_counted();
    std::cout << "rhs calls (tol = 1e-6) = " << counted.count() << "\n";

    return 0;
}
