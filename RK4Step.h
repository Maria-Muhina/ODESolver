#pragma once

#include "OdeSystem.h"
#include "OdeTypes.h"

namespace ode {

    State rk4_step(
        const OdeSystem& system,
        double t,
        const State& x,
        double h
    );

}