#include <stdexcept>
#include <cmath>
#include "OdeValidation.h"

namespace ode {
    void validate_time_parameters(double t0, double t_end, double h) {
        if (!std::isfinite(h) || !std::isfinite(t0) || !std::isfinite(t_end)) {
            throw std::invalid_argument("All parameters must be finite");
        }
        if (h <= 0) {
            throw std::invalid_argument("Time step must be positive");
        }
        if (t0 > t_end) {
            throw std::invalid_argument("Start time cannot be greater than end time");
        }
    }
}