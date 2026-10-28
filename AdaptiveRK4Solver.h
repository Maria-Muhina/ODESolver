#pragma once

#include "OdeSystem.h"
#include "OdeTypes.h"

namespace ode {

    class AdaptiveRK4Solver{
        public:
            Trajectory solve(
                const OdeSystem& system,
                double t0,
                const State& x0,
                double t_end,
                double h_initial,
                double tol
            ) const;
    };
}