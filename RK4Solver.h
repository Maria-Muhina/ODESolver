#pragma once

#include "OdeTypes.h"
#include "OdeSystem.h"

namespace ode
{
    class RK4Solver
    {
    public:
        Trajectory solve(
            const OdeSystem& system,
            double t0,
            const State& x0,
            double t_end,
            double h
        ) const;
};
    
} 
