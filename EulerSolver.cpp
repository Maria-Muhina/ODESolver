#include <algorithm>
#include <stdexcept>
#include "EulerSolver.h"
#include "OdeValidation.h"

ode::Trajectory ode::EulerSolver::solve(
    const OdeSystem& system, double t0, const ode::State& x0, double t_end, double h) const {
        ode::validate_time_parameters(t0, t_end, h);
        ode::Trajectory result;
        result.push_back({t0, x0});
        double t = t0;
        ode::State x = x0;
        ode::State dxdt(x0.size());
        std::size_t n = 0;
        
        while (t < t_end) {
            n++;
            double candidate_t = t0 + static_cast<double>(n) * h;
            double next_t = std::min(candidate_t, t_end);
            if (next_t <= t) {
                throw std::runtime_error("Time step does not advance time");
            }
            double step = next_t - t;
            system.rhs(t, x, dxdt);
            x = x + step * dxdt;
            t = next_t;
            result.push_back({t, x});
        }

        return result;
    }