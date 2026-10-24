#include <algorithm>
#include <stdexcept>
#include "RK4Solver.h"
#include "OdeValidation.h"

ode::Trajectory ode::RK4Solver::solve(
    const ode::OdeSystem& system, double t0, const ode::State& x0, double t_end, double h) const {
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
            ode::State k1 = dxdt; 
            system.rhs(t + step / 2, x + step / 2 * k1, dxdt);
            ode::State k2 = dxdt;
            system.rhs(t + step / 2, x + step / 2 * k2, dxdt);
            ode::State k3 = dxdt;
            system.rhs(t + step, x + step * k3, dxdt);
            ode::State k4 = dxdt;

            x += step / 6 * (k1 + 2*k2 + 2*k3 + k4);
            t = next_t;

            result.push_back({t, x});
        }
        return result;
    }
