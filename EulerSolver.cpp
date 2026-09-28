#include "EulerSolver.h"

std::vector<ode::SolutionPoint> ode::EulerSolver::solve(
    const OdeSystem& system, double t0, double x0, double t_end, double h) const {
            std::vector<SolutionPoint> result;
            result.push_back({t0, x0});
            double t = t0;
            double x = x0;
            while (t < t_end) {
                x = x + h * system.rhs(t, x);
                t = t + h;
                result.push_back({t, x});
            }

            return result;
        }