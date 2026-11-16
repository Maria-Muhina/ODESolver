#pragma once

#include "OdeTypes.h"
#include "OdeSystem.h"

namespace ode {
    void validate_time_parameters(
        double t0,
        double t_end,
        double h
    );

    void checked_rhs(
    const OdeSystem& system,
    double t,
    const State& x,
    State& dxdt
    );
}
