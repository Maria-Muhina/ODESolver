#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>

Eigen::VectorXd rhs(double t, const Eigen::VectorXd& y) {
    assert(y.size() == 2);

    Eigen::VectorXd result(2);
    result << y[1], -1 * y[0];

    (void)t;
    return result;
}


int main() {
    Eigen::VectorXd x(2);
    x << 1.0, 2.0;
    
    double norm = x.norm();

    std::cout << "norm = " << norm << "\n";

    std::cout << "x size = " << x.size() << "\n";

    Eigen::VectorXd y(2);
    y << 3.0, 2.0;

    Eigen::VectorXd dy = rhs(0.0, y);

    std::cout << "y = " << y << "\n";

    std::cout << "dy = " << dy << "\n";

    Eigen::VectorXd example_1(2);
    example_1 << 3, 2;
    auto result_1 = rhs(0.0, example_1);
    std::cout << "result_1 = " << result_1 << "\n";

    assert(std::abs(result_1[0] - 2) < 1e-14);
    assert(std::abs(result_1[1] - -3) < 1e-14);

    Eigen::VectorXd example_2(2);
    example_2 << -1, 4;
    auto result_2 = rhs(0.5, example_2);
    std::cout << "result_2 = " << result_2 << "\n";

    assert(std::abs(result_2[0] - 4) < 1e-14);
    assert(std::abs(result_2[1] - 1) < 1e-14);

    Eigen::VectorXd example_3(2);
    example_3 << 0, 0;
    auto result_3 = rhs(1.0, example_3);
    std::cout << "result_3 = " << result_3 << "\n";

    assert(std::abs(result_3[0] - 0) < 1e-14);
    assert(std::abs(result_3[1] - 0) < 1e-14);
}
