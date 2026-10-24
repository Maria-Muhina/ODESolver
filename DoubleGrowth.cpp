#include "DoubleGrowth.h"

void ode::DoubleGrowth::rhs(
            double t,
            const State& x,
            State& dxdt
        ) const {
    dxdt(0) = 2*x(0);
};