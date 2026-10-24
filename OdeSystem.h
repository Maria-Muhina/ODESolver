#pragma once

#include "OdeTypes.h"

namespace ode {
class OdeSystem {
    public:
        virtual void rhs(
            double t,
            const State& x,
            State& dxdt
        ) const = 0;    

        virtual ~OdeSystem() = default;
};

void check_system(const OdeSystem& system, double t, const State& x, const State& expected);

}