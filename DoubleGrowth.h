#pragma once

#include "OdeSystem.h"

namespace ode {
    class DoubleGrowth : public OdeSystem {
        public:
            double rhs(double t, double x) const override;
    };
}