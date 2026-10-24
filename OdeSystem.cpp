#include <cassert>
#include <cmath>
#include "OdeSystem.h"

namespace ode
{
    void check_system(const OdeSystem &system, double t, const State &x, const State &expected)
    {
        assert(x.size() == expected.size());
        State dxdt(x.size());
        system.rhs(t, x, dxdt);
        for (Eigen::Index i = 0; i < x.size(); i++)
        {
            assert(std::abs(dxdt(i) - expected(i)) < 1e-14);
        }
    };

}
