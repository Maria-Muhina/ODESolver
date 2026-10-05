#include "RK4Solver.h"

std::vector<ode::SolutionPoint> ode::RK4Solver::solve(
    const ode::OdeSystem& system, double t0, double x0, double t_end, double h) const {
        std::vector<ode::SolutionPoint> result;
        result.push_back({t0, x0});
        double t = t0;
        double x = x0;
        while (t < t_end) {
            double step = h;
            if (t + h > t_end) {
                step = t_end - t;
            }

            double k1 = system.rhs(t, x);
            double k2 = system.rhs(t + step / 2, x + step / 2 * k1);
            double k3 = system.rhs(t + step / 2, x + step / 2 * k2);
            double k4 = system.rhs(t + step, x + step * k3);

            x += step / 6 * (k1 + 2*k2 + 2*k3 + k4);
            t += step;
            result.push_back({t, x});
        }
        return result;
    }
