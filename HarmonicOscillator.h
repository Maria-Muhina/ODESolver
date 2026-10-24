#pragma once

#include "OdeSystem.h"

namespace ode {
class HarmonicOscillator : public OdeSystem {
    public:
        void rhs(
            double t,
            const State& x,
            State& dxdt
        ) const override;
    };

}