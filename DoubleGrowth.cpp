#include "DoubleGrowth.h"

double ode::DoubleGrowth::rhs(double t, double x) const {
    return 2*x;
};