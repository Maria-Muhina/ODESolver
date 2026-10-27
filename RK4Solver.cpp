#include <algorithm>
#include <stdexcept>
#include "RK4Solver.h"
#include "OdeValidation.h"
#include "RK4Step.h"

ode::Trajectory ode::RK4Solver::solve(
    const ode::OdeSystem& system, double t0, const ode::State& x0, double t_end, double h) const {
        ode::validate_time_parameters(t0, t_end, h);
        ode::Trajectory result;
        result.push_back({t0, x0});
        double t = t0;
        ode::State x = x0;
        std::size_t n = 0;

        while (t < t_end) {
            n++;
            double candidate_t = t0 + static_cast<double>(n) * h;
            double next_t = std::min(candidate_t, t_end);
            if (next_t <= t) {
                throw std::runtime_error("Time step does not advance time");
            }
            double step = next_t - t;

            x = ode::rk4_step(system, t, x, step);
            t = next_t;

            result.push_back({t, x});
        }
        return result;
    }
