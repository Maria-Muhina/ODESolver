#pragma once

#include <Eigen/Core>
#include <vector>

namespace ode
{
    using State = Eigen::VectorXd;

    struct SolutionPoint {
        double t;
        State x;
    };

    using Trajectory = std::vector<SolutionPoint>;
} // namespace ode
