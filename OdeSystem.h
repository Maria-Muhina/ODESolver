#pragma once

namespace ode {
    class OdeSystem {
        public:
            virtual double rhs(double t, double x) const = 0;
    };

    void check_system(const OdeSystem& system, double t, double x, double expected);
}