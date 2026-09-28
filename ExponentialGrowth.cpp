#include "ExponentialGrowth.h"

double ode::ExponentialGrowth::rhs(double t, double x) const {
    return x;
}