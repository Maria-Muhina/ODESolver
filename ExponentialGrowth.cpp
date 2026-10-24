#include "ExponentialGrowth.h"

void ode::ExponentialGrowth::rhs(
            double t,
            const State& x,
            State& dxdt
        ) const {
            dxdt(0) = x(0);
        }
