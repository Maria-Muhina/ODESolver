#include "HarmonicOscillator.h"

void ode::HarmonicOscillator::rhs(
            double t,
            const State& x,
            State& dxdt
        ) const {
    dxdt(0) = x(1);
    dxdt(1) = -x(0);
}