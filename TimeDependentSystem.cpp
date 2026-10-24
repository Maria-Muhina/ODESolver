#include "TimeDependentSystem.h"

void ode::TimeDependentSystem::rhs(
            double t,
            const State& x,
            State& dxdt
        ) const {
    dxdt(0) = t;
}