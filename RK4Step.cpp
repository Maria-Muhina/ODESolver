#include "RK4Step.h"

ode::State ode::rk4_step(const ode::OdeSystem& system, double t, const ode::State& x, double h) {
    ode::State dxdt(x.size());

    system.rhs(t, x, dxdt);
    ode::State k1 = dxdt; 
    system.rhs(t + h / 2, x + h / 2 * k1, dxdt);
    ode::State k2 = dxdt;
    system.rhs(t + h / 2, x + h / 2 * k2, dxdt);
    ode::State k3 = dxdt;
    system.rhs(t + h, x + h * k3, dxdt);
    ode::State k4 = dxdt;

    return x + h / 6.0 * (k1 + 2*k2 + 2*k3 + k4);
}
