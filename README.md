# Numerical Dynamics in C++

A small C++ project for implementing and testing numerical methods for ordinary differential equations (ODEs).

The project implements three explicit ODE solvers: Forward Euler, classical fourth-order Runge–Kutta (RK4) with a fixed nominal step, and adaptive RK4 with local error estimation by step doubling.

The implementation is intentionally small. The goal at this stage is to establish a clean numerical-solver foundation before adding more advanced numerical methods and engineering infrastructure.

---

## Overview

The project considers vector-valued initial-value problems of the form

```math
x'(t) = f(t,x), \qquad x(t_0) = x_0, \qquad x(t) \in \mathbb{R}^n.
```

When $n=1$, the problem reduces to a scalar ODE. A numerical solver takes an ODE system, an initial condition, a final time, and an initial or nominal step size, and returns the computed solution as a trajectory of time/state pairs.

The main components are:

- `OdeSystem` — interface for the right-hand side of an ODE;
- `EulerSolver` — Forward Euler implementation;
- `RK4Solver` — classical fourth-order Runge–Kutta implementation with a fixed nominal step;
- `RK4Step` — reusable implementation of a single classical RK4 step;
- `AdaptiveRK4Solver` — adaptive RK4 implementation using step doubling;
- concrete ODE systems used for testing;
- GoogleTest-based automated verification.

The fixed-step solvers share the same basic interface. The adaptive solver adds an absolute local error tolerance and interprets the supplied step size as the initial step.

---

## Quick start

The following example integrates the scalar test problem $x'(t)=x(t)$, $x(0)=1$, from $t=0$ to $t=1$ using adaptive RK4:

```cpp
#include "AdaptiveRK4Solver.h"
#include "ExponentialGrowth.h"

int main() {
    ode::ExponentialGrowth system;
    ode::AdaptiveRK4Solver solver;

    ode::State x0(1);
    x0 << 1.0;

    ode::Trajectory solution = solver.solve(
        system,
        0.0,    // initial time
        x0,     // initial state
        1.0,    // final time
        0.2,    // initial step size
        1e-6    // absolute local error tolerance
    );

    double final_value = solution.back().x(0);
    (void)final_value; // use the computed value as needed
}
```

`solution` contains the initial state and all accepted integration points, each with a time `t` and a state vector `x`. This is a usage example, not a separate executable target configured in the supplied `CMakeLists.txt`.

---

## Project structure

The project currently uses a flat source layout at the repository root, rather than separate `src/`, `include/`, and `tests/` directories. The following is a grouped inventory of the relevant files (not a directory hierarchy):

```text
Core numerical API:
  OdeTypes.h
  OdeSystem.h / OdeSystem.cpp
  OdeValidation.h / OdeValidation.cpp
  EulerSolver.h / EulerSolver.cpp
  RK4Step.h / RK4Step.cpp
  RK4Solver.h / RK4Solver.cpp
  AdaptiveRK4Solver.h / AdaptiveRK4Solver.cpp

Test systems:
  ExponentialGrowth.h / ExponentialGrowth.cpp
  DoubleGrowth.h / DoubleGrowth.cpp
  TimeDependentSystem.h / TimeDependentSystem.cpp
  HarmonicOscillator.h / HarmonicOscillator.cpp

Automated tests:
  EulerSolverTests.cpp
  RK4SolverTests.cpp
  AdaptiveRK4SolverTests.cpp

Project and experiments:
  CMakeLists.txt
  README.md
  main.cpp
  Eigen_test.cpp
  Counter.h / Counter.cpp
  .gitignore
```

`Counter.*` is from the initial skeleton and is not part of the numerical solver targets.

| Component | Responsibility |
|---|---|
| `OdeTypes` | Defines `State`, `SolutionPoint`, and `Trajectory` |
| `OdeSystem` | Defines the interface for evaluating the ODE right-hand side |
| `OdeValidation` | Validates integration time parameters |
| `EulerSolver` | Integrates the system using Forward Euler |
| `RK4Solver` | Integrates the system using classical RK4 with a fixed nominal step |
| `RK4Step` | Computes one RK4 step and returns a new state |
| `AdaptiveRK4Solver` | Selects steps using RK4 step doubling and an absolute error tolerance |
| Concrete ODE systems | Provide test problems for verifying numerical methods |

`State` is represented by `Eigen::VectorXd`. The `OdeSystem::rhs(double t, const State& x, State& dxdt)` interface writes the derivative into a caller-provided output vector. The Forward Euler solver reuses a derivative buffer during integration. The shared `rk4_step()` function allocates a local derivative buffer for each call. The output vector must have the dimension expected by the system.

`SolutionPoint` stores a time and its corresponding state; `Trajectory` is a sequence of solution points. All three solvers return a trajectory; the adaptive solver stores accepted points only.

`main.cpp` provides manual checks and examples; `Eigen_test.cpp` is an Eigen experimentation target. Automated regression tests cover Euler, fixed-step RK4, the shared RK4 step, and adaptive RK4. The one-step RK4 test is part of the RK4 test source, not a separate test source file.

Generated build and test artifacts are not part of the project source code.

---

## Mathematical foundation

### Problem formulation

The solver works with an initial-value problem

```math
x'(t) = f(t,x), \qquad x(t_0) = x_0,
```

where $x(t), x_0 \in \mathbb{R}^n$ and $f(t,x) \in \mathbb{R}^n$.

Here:

- $t$ is the independent variable;
- $x(t)$ is the state vector;
- $f(t,x)$ is the right-hand side of the differential equation;
- $t_0$ and $x_0$ define the initial condition.

The numerical solution is computed at a sequence of time points

```math
t_0,\ t_1,\ \ldots,\ t_N.
```

For fixed-step integration, the nominal step size is $h$. For adaptive RK4, $h$ is the initial step size and changes according to the estimated local error. In all solvers, a proposed step is shortened when necessary to avoid passing `t_end`.

### Primary test problem

The main test problem used throughout the project is

```math
x'(t)=x(t), \qquad x(0)=1.
```

Its exact solution is

```math
x(t)=e^t.
```

This problem is useful because the exact solution is known, allowing the numerical solution to be compared directly with the analytical one. It is also used for convergence experiments for both Euler and RK4.

### Additional test systems

The project also contains several other ODE systems:

- $x'(t)=2x(t)$, used to verify that the solvers work with a different growth rate.
- $x'(t)=t$, $x(0)=1$, whose exact solution is

  ```math
  x(t)=1+\frac{t^2}{2}.
  ```

  This problem verifies correct handling of time-dependent right-hand sides.
- A two-dimensional harmonic oscillator:

  ```math
  \begin{cases}
  q'(t)=v(t),\\
  v'(t)=-q(t),
  \end{cases}
  ```

  with initial conditions

  ```math
  q(0)=1, \qquad v(0)=0.
  ```

  The exact solution is

  ```math
  q(t)=\cos(t), \qquad v(t)=-\sin(t).
  ```

  This system is used to verify that both fixed-step and adaptive RK4 handle vector-valued states and coupled ODEs.

Using several systems helps avoid testing the solver against only one specific form of $f(t,x)$.

---

## Numerical methods

### Forward Euler

The Forward Euler method evaluates the derivative at the beginning of each step:

```math
x_{n+1}=x_n+h f(t_n,x_n).
```

Here $h$ denotes the actual length of the step; the final step may be shorter than the nominal step size.

For the primary test problem $x'=x$, this becomes

```math
x_{n+1}=x_n+h x_n.
```

Under standard regularity assumptions, the method has first-order global convergence:

```math
E(h)\approx Ch.
```

### Classical RK4

The classical fourth-order Runge–Kutta method evaluates the right-hand side at four stages:

```math
\begin{aligned}
k_1 &= f(t_n,x_n),\\
k_2 &= f\left(t_n+\frac h2,\ x_n+\frac h2 k_1\right),\\
k_3 &= f\left(t_n+\frac h2,\ x_n+\frac h2 k_2\right),\\
k_4 &= f(t_n+h,\ x_n+h k_3).
\end{aligned}
```

The solution is then updated using

```math
x_{n+1}=x_n+\frac h6\left(k_1+2k_2+2k_3+k_4\right).
```

Here, too, $h$ is the actual length of the step. Under standard regularity assumptions and for sufficiently small steps, RK4 has fourth-order global convergence:

```math
E(h)\approx Ch^4.
```

Fixed-step RK4 uses the same solver interface as Euler, making the methods directly comparable on the same ODE systems. The four-stage calculation is implemented once in `rk4_step()` and reused by both RK4 solvers.

---

### Adaptive RK4 with step doubling

`AdaptiveRK4Solver` uses the classical RK4 step formula, but adapts the step size based on an estimate of the local error.

For each attempted step of length $h$, the solver computes:

- `x_full`: one RK4 step of length $h$;
- `x_half1`: one RK4 step of length $h/2$;
- `x_half2`: another RK4 step of length $h/2$, starting from `x_half1` at time $t+h/2$.

For sufficiently small steps, the estimate of the local error in the two-half-step approximation is

```math
e = \frac{\|x_{\mathrm{half2}}-x_{\mathrm{full}}\|_2}{15}.
```

The attempted step is accepted when $e \leq \mathrm{tol}$, where `tol` is a positive absolute local error tolerance. On acceptance, the solver advances time by the attempted step and stores `x_half2`. On rejection, the time and state remain unchanged and the step is retried with a smaller size. Rejected attempts are not stored in the trajectory.

The next proposed step size is calculated from the *actual attempted* step length:

```math
h_{\mathrm{new}} = s\,h_{\mathrm{used}}
\left(\frac{\mathrm{tol}}{e}\right)^{1/5},
\qquad s=0.9.
```

The multiplicative step-size factor is clamped to the range $[0.2, 5.0]$. When the error estimate is zero, the factor is set to $5.0$ instead of dividing by zero. The actual attempted step is also limited by the remaining integration interval.

The exponent $1/5$ follows from the leading-order local truncation error of classical RK4. The error control is based on an **absolute local estimate**, not a guarantee on the accumulated global error.

Each attempted adaptive step performs three RK4 steps, or **12 evaluations of `rhs()`**. A fixed RK4 step performs four. Rejected adaptive attempts also incur this cost.

---

## Verification and testing

Automated tests cover method-specific behavior as well as the reusable RK4 step and the adaptive controller.

The test suite currently contains **24 GoogleTest test cases** (9 Euler, 8 fixed-step RK4, 1 reusable RK4 step, and 6 adaptive RK4). All 24 passed in the latest user-reported run.

The tests cover:

- integration of the primary test problem;
- reaching the requested final time, including cases where floating-point rounding would otherwise add an extra step;
- different initial conditions and right-hand sides;
- time-dependent right-hand sides;
- numerical accuracy against known exact solutions;
- decreasing error with a smaller step size and observed convergence order;
- integration of a two-dimensional harmonic oscillator;
- rejection of invalid integration parameters;
- correct handling of a zero-length integration interval;
- adaptive rejection of a large initial step and increasing time at accepted points;
- zero-error step-size growth and a vector-valued adaptive oscillator test;
- invalid adaptive tolerances and inability to advance floating-point time.

### Convergence experiments

For a method with global error

```math
E(h)\approx Ch^p,
```

reducing the step size by a factor of two gives approximately

```math
\frac{E(h)}{E(h/2)}\approx 2^p.
```

The experimental order can therefore be estimated as

```math
p\approx\frac{\log(E(h)/E(h/2))}{\log 2}.
```

#### Euler

For the primary test problem $x'(t)=x(t)$, $x(0)=1$, the endpoint error was evaluated at $t=0.3$ for these nominal step sizes:

```text
h = 0.1
h = 0.05
h = 0.025
```

The measured absolute errors were approximately:

```text
h       E(h)
0.1     1.88588e-2
0.05    9.76317e-3
0.025   4.96998e-3
```

The corresponding experimental orders were approximately:

```text
0.1  → 0.05   : p ≈ 0.95
0.05 → 0.025  : p ≈ 0.97
```

The estimates approach the expected first-order behavior as the step size decreases.

#### RK4

For the same primary test problem, the endpoint error was evaluated at $t=0.3$ for:

```text
h = 0.2
h = 0.1
h = 0.05
h = 0.025
```

The measured absolute errors were approximately:

```text
h       E(h)
0.2     3.15174e-6
0.1     3.10513e-7
0.05    2.02313e-8
0.025   1.29105e-9
```

The corresponding experimental orders were approximately:

```text
0.2  → 0.1    : p ≈ 3.34
0.1  → 0.05   : p ≈ 3.94
0.05 → 0.025  : p ≈ 3.97
```

The first estimate is lower because the coarser step size is not yet fully in the asymptotic convergence regime. As the step size decreases, the measured order approaches the expected value of 4.

For $h=0.2$, the interval $[0,0.3]$ is not covered by an integer number of full steps, so the implementation uses a shortened final step. This is part of the solver's normal final-step handling.

The observed increase of these particular order estimates is a property of this experiment; it is not assumed to be a general requirement for every RK4 convergence experiment.

---

### Adaptive RK4 accuracy and evaluation cost

For $x'(t)=x(t)$ with $x(0)=1$, integrated on $[0,1]$ with initial step $0.2$, the observed adaptive results were:

| Absolute tolerance | Accepted steps | Absolute endpoint error |
|---:|---:|---:|
| `1e-3` | 2 | `3.33360e-4` |
| `1e-6` | 5 | `3.85803e-6` |
| `1e-9` | 17 | `1.74269e-8` |

Tighter local tolerances increased the number of accepted steps and reduced the measured endpoint error in this experiment. The endpoint error may exceed `tol`, because `tol` controls an estimate of **local** rather than global error.

A separate exploratory comparison measured right-hand-side evaluations for the same equation and final time:

| Method | Step / tolerance | `rhs()` evaluations | Absolute endpoint error |
|---|---|---:|---:|
| Fixed RK4 | `h = 0.2` | 20 | `3.06919e-5` |
| Fixed RK4 | `h = 0.1` | 40 | `2.08432e-6` |
| Fixed RK4 | `h = 0.05` | 80 | `1.35803e-7` |
| Adaptive RK4 | `h_initial = 0.2`, `tol = 1e-6` | 60 | `3.85803e-6` |

In this short, smooth problem, fixed-step RK4 at $h=0.1$ achieved a smaller endpoint error with fewer `rhs()` evaluations than the adaptive solver. This does not imply that fixed steps are generally more efficient: step doubling has extra computational overhead, while adaptive stepping can be useful when solution behavior changes substantially over an interval.

The evaluation counts come from a temporary counting experiment, not a permanent regression test or a full runtime benchmark.

---

## Building

### Requirements

- CMake **3.16 or newer**.
- A compatible C++ compiler and build tool (the project enables `/W4` for MSVC or `-Wall -Wextra -Wpedantic` for other compilers).
- Network access on the first configuration so CMake `FetchContent` can download **Eigen 3.4.0** and **GoogleTest 1.14.0**.

**C++ standard:** the supplied `CMakeLists.txt` does not explicitly set `CMAKE_CXX_STANDARD` or `target_compile_features`. The exact minimum supported standard has therefore not been established by the project's build configuration. Pinning and testing a standard is a future reproducibility improvement.

### Configure and build

Run these commands from the repository root. First, configure the project:

```bash
cmake -S . -B build
```

Then build:

```bash
cmake --build build
```

The CMake configuration defines three executable targets: `main`, `ODESolverTests`, and `Eigen_test`.

---

## Running tests

Run the project's ODE tests (excluding tests registered by the Eigen dependency):

```bash
ctest --test-dir build -R '^(EulerSolverTest|RK4SolverTest|RK4StepTest|AdaptiveRK4SolverTest)\.' --output-on-failure
```

The test executable can also be run directly:

```bash
./build/ODESolverTests
```

To run a specific GoogleTest test:

```bash
./build/ODESolverTests --gtest_filter=RK4SolverTest.SolvesDoubleGrowth
```

---

## Current limitations

The current implementation deliberately has a limited scope. It supports:

- scalar and vector-valued ODE systems through `Eigen::VectorXd`;
- fixed nominal step sizes with explicit Euler and classical RK4;
- adaptive classical RK4 using step doubling and an absolute local error tolerance;
- forward integration with a shortened final step;
- automated verification with GoogleTest.

Current limitations include:

- Only forward integration is supported (`t_end >= t0`, `h > 0`).
- The entire computed trajectory is stored in memory.
- Adaptive RK4 currently uses only a scalar absolute local tolerance; relative and component-wise scaled tolerances are not implemented.
- Stiff ODE methods, event detection, and backward integration are not implemented.
- State dimension consistency and state-value finiteness are not fully validated.
- Time-step advancement is checked and a runtime error is raised if a proposed step cannot advance the floating-point time. Extreme floating-point scales, very small tolerances, and overflow or underflow in step-size control are not comprehensively handled.
- A production-level safeguard against excessive rejected attempts and comprehensive validation of error-controller intermediates are not implemented.
- Systematic benchmarking, profiling, continuous integration, and extended public API documentation are not yet implemented.

These features are intentionally outside the current implementation scope.

---

## Development direction

The project is being developed incrementally, with emphasis on:

- numerical correctness;
- reproducible verification;
- clear solver interfaces;
- maintainable C++ code;
- automated testing;
- progressively stronger software-engineering practices.

The Euler, fixed-step RK4, and adaptive RK4 implementations form the numerical foundation for later work on solver robustness, accuracy controls, and scientific-software engineering.
