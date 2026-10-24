#pragma once

#include "OdeSystem.h"

namespace ode {
class DoubleGrowth : public OdeSystem {
    public:
        virtual void rhs(
            double t,
            const State& x,
            State& dxdt
        ) const override;
    };

}