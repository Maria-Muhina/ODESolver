#include <cassert>
#include "OdeSystem.h"

namespace ode {
    void check_system(const ode::OdeSystem& system, double t, double x, double expected) {
        assert(system.rhs(t, x) == expected);
    };
}
