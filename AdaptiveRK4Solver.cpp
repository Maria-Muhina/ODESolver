#include <stdexcept>
#include <cmath>
#include <algorithm>
#include "AdaptiveRK4Solver.h"
#include "OdeValidation.h"
#include "RK4Step.h"

ode::Trajectory ode::AdaptiveRK4Solver::solve(
    const ode::OdeSystem& system, double t0, const ode::State& x0, 
        double t_end, double h_initial, double tol) const {

        ode::validate_time_parameters(t0, t_end, h_initial);

        if (!std::isfinite(tol)) {
            throw std::invalid_argument("Tolerance must be finite");
        }
        if (tol <= 0.0) {
            throw std::invalid_argument("Tolerance must be positive");
        }

        ode::Trajectory result;
        result.push_back({t0, x0});

        double t = t0;
        ode::State x = x0;
        double h = h_initial;

        while (t < t_end) {
            double step = std::min(h, t_end - t);
            if (t + step <= t) {
                throw std::runtime_error("Time step does not advance time");
            }

            ode::State x_full = rk4_step(system, t, x, step);
            ode::State x_half1 = rk4_step(system, t, x, step / 2);
            ode::State x_half2 = rk4_step(system, t + step / 2, x_half1, step / 2);

            double e = (x_half2 - x_full).norm() / 15.0;

            
            if (e <= tol) {
                x = x_half2;
                t = t + step;
                result.push_back({t, x});
            }

            double factor;

            if (!std::isfinite(e)) {
                throw std::runtime_error("Non-finite error estimate");
            }

            if (e > 0.0) {
                factor = 0.9 * std::pow(tol / e, 1.0 / 5.0);
                factor = std::clamp(factor, 0.2, 5.0);
            } else {
                factor = 5.0;
            }

            h = factor * step;

        }


        return result;
    }
