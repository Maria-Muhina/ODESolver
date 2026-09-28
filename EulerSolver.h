#pragma once

#include <vector>
#include "OdeSystem.h"


namespace ode {

struct SolutionPoint {
    double t;
    double x;
};

class EulerSolver {
    public:
        std::vector<SolutionPoint> solve(
            const OdeSystem& system,
            double t0,
            double x0,
            double t_end,
            double h
        ) const;
};

}