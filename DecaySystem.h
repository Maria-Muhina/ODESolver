#pragma once

#include "OdeSystem.h"

namespace ode {

class DecaySystem : public OdeSystem {
public:
    void rhs(
        double /*t*/,
        const State& x,
        State& dxdt
    ) const override {
        dxdt(0) = -10.0 * x(0);
    }
};

} // namespace ode