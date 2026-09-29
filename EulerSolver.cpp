#include "EulerSolver.h"

std::vector<ode::SolutionPoint> ode::EulerSolver::solve(
    const OdeSystem& system, double t0, double x0, double t_end, double h) const {
            std::vector<SolutionPoint> result;
            result.push_back({t0, x0});
            double t = t0;
            double x = x0;
            while (t < t_end) {
                double step = h;
                if (t + h <= t_end) {
                    step = h;
                }else {
                    step = t_end - t;
                }
                x = x + step * system.rhs(t, x);
                t = t + step;
                result.push_back({t, x});
            }

            return result;
        }