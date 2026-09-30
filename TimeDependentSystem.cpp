#include "TimeDependentSystem.h"

double ode::TimeDependentSystem::rhs(double t, double x) const {
    return t;
}