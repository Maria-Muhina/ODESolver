#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

Eigen::VectorXd rhs(double t, const Eigen::VectorXd& y) {
    assert(y.size() == 2);

    Eigen::VectorXd result(2);
    result << y[1], -1 * y[0];

    (void)t;
    return result;
}

Eigen::VectorXd euler_step(
    double t,
    const Eigen::VectorXd& y,
    double h
) {
    Eigen::VectorXd result(y.size());
    result = y + h * rhs(t, y);

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


    Eigen::VectorXd y0(2);
    y0 << 3.0, 2.0;

    auto result_euler = euler_step(0.0, y0, 0.1);
    std::cout << "result_euler = " << result_euler << "\n";

    assert(std::abs(result_euler[0] - 3.2) < 1e-14);
    assert(std::abs(result_euler[1] - 1.7) < 1e-14);
    assert(std::abs(y0[0] - 3) < 1e-14);
    assert(std::abs(y0[1] - 2) < 1e-14);


    auto result_euler_2 = euler_step(0.0, y0, 0.0);
    std::cout << "result_euler_2 = " << result_euler_2 << "\n";

    assert(std::abs(result_euler_2[0] - 3) < 1e-14);
    assert(std::abs(result_euler_2[1] - 2) < 1e-14);

    Eigen::VectorXd y0_0(2);
    y0_0 << 0.0, 0.0;

    auto result_euler_3 = euler_step(0.0, y0_0, 0.1);
    std::cout << "result_euler_3 = " << result_euler_3 << "\n";

    assert(std::abs(result_euler_3[0] - 0.0) < 1e-14);
    assert(std::abs(result_euler_3[1] - 0.0) < 1e-14);

    Eigen::VectorXd y_result(2);
    double t = 0.0;
    double h = 0.1;
    y_result = y;
    for (int i = 0; i < 10; i ++){
        y_result = euler_step(t, y_result, h);
        t += h;
    }
    std::cout << "t = " << t << "\n";
    std::cout << "y_result = " << y_result << "\n";


    std::vector<double> step_sizes = {0.1, 0.05, 0.025};
    double T = 1.0;
    std::vector<double> errors(3);
    int k = 0;
    for (double h : step_sizes) {
        Eigen::VectorXd y0_loop(2);
        y0_loop << 3, 2;
        Eigen::VectorXd y_result_loop = y0_loop;
        double t = 0;
        int N = static_cast<int>(T / h);
        for (int i = 0; i < N; i ++) {
            y_result_loop = euler_step(t, y_result_loop, h);
            t += h;
        }
        double q_err = y_result_loop[0] - (3*std::cos(1) + 2*std::sin(1));
        double v_err = y_result_loop[1] - (-3*std::sin(1) + 2*std::cos(1));
        double err = std::sqrt(q_err*q_err + v_err*v_err);
        std::cout << "err_" << h << " = " << err << "\n";
        errors[k] = err;
        k += 1;
    }
    
    double p_1 = std::log(errors[0]/errors[1]) / std::log(2);
    double p_2 = std::log(errors[1]/errors[2]) / std::log(2);
    std::cout << "p_1 = " << p_1 << "\n";
    std::cout << "p_2 = " << p_2 << "\n";
}
