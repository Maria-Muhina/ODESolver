# Numerical Dynamics in C++

A small C++ project for implementing and testing numerical methods for ordinary differential equations (ODEs).

The project currently focuses on two explicit time-stepping methods: the Forward Euler method and the classical fourth-order Runge–Kutta method (RK4).

The implementation is intentionally small. The goal at this stage is to establish a clean numerical-solver foundation before adding more advanced methods and engineering infrastructure.

---

## Overview

The project considers vector-valued initial-value problems of the form

$$
x'(t) = f(t,x), \qquad x(t_0) = x_0, \qquad x(t) \in \mathbb{R}^n.
$$

When $n=1$, the problem reduces to a scalar ODE. A numerical solver takes an ODE system, an initial condition, a final time, and a nominal step size, and returns the computed solution at the integration points.

The main components are:

- `OdeSystem` — interface for the right-hand side of an ODE;
- `EulerSolver` — Forward Euler implementation;
- `RK4Solver` — classical fourth-order Runge–Kutta implementation;
- concrete ODE systems used for testing;
- GoogleTest-based automated verification.

Both solvers use the same basic interface so that different numerical methods can be applied to the same initial-value problem.

---

## Project structure

The main numerical components are organized as follows:

```text
.
├── OdeTypes.h
├── OdeSystem.h / .cpp
├── OdeValidation.h / .cpp
├── EulerSolver.h / .cpp
├── RK4Solver.h / .cpp
├── ExponentialGrowth.h / .cpp
├── DoubleGrowth.h / .cpp
├── TimeDependentSystem.h / .cpp
├── HarmonicOscillator.h / .cpp
├── EulerSolverTests.cpp
├── RK4SolverTests.cpp
├── main.cpp
├── CMakeLists.txt
└── .gitignore
```

| Component | Responsibility |
|---|---|
| `OdeTypes` | Defines `State`, `SolutionPoint`, and `Trajectory` |
| `OdeSystem` | Defines the interface for evaluating the ODE right-hand side |
| `OdeValidation` | Validates integration time parameters |
| `EulerSolver` | Integrates the system using Forward Euler |
| `RK4Solver` | Integrates the system using classical RK4 |
| Concrete ODE systems | Provide test problems for verifying numerical methods |

`State` is represented by `Eigen::VectorXd`. The `OdeSystem::rhs(double t, const State& x, State& dxdt)` interface writes the derivative into a caller-provided output vector. Solvers allocate the derivative buffer and reuse it during integration. The output vector must have the dimension expected by the system.

`SolutionPoint` stores a time and its corresponding state; `Trajectory` is a sequence of solution points. Both solvers return the complete computed trajectory.

`main.cpp` provides manual checks and examples. The automated regression tests are in `EulerSolverTests.cpp` and `RK4SolverTests.cpp`. Other files from the initial project skeleton, such as `Counter.*`, are not part of the numerical-solver design.

Generated build and test artifacts are not part of the project source code.

---

## Mathematical foundation

### Problem formulation

The solver works with an initial-value problem

$$
x'(t) = f(t,x), \qquad x(t_0) = x_0,
$$

where $x(t), x_0 \in \mathbb{R}^n$ and $f(t,x) \in \mathbb{R}^n$.

Here:

- $t$ is the independent variable;
- $x(t)$ is the state vector;
- $f(t,x)$ is the right-hand side of the differential equation;
- $t_0$ and $x_0$ define the initial condition.

The numerical solution is computed at a sequence of time points

$$
t_0,\ t_1,\ \ldots,\ t_N.
$$

The nominal step size is $h$. If the requested final time is not reached exactly by a full step, the final step is shortened to end at `t_end`.

### Primary test problem

The main test problem used throughout the project is

$$
x'(t)=x(t), \qquad x(0)=1.
$$

Its exact solution is

$$
x(t)=e^t.
$$

This problem is useful because the exact solution is known, allowing the numerical solution to be compared directly with the analytical one. It is also used for convergence experiments for both Euler and RK4.

### Additional test systems

The project also contains several other ODE systems:

- $x'(t)=2x(t)$, used to verify that the solvers work with a different growth rate.
- $x'(t)=t$, $x(0)=1$, whose exact solution is

  $$
  x(t)=1+\frac{t^2}{2}.
  $$

  This problem verifies correct handling of time-dependent right-hand sides.
- A two-dimensional harmonic oscillator:

  $$
  \begin{cases}
  q'(t)=v(t),\\
  v'(t)=-q(t),
  \end{cases}
  $$

  with initial conditions

  $$
  q(0)=1, \qquad v(0)=0.
  $$

  The exact solution is

  $$
  q(t)=\cos(t), \qquad v(t)=-\sin(t).
  $$

  This system is used to verify that RK4 handles vector-valued states and integrates coupled ODEs.

Using several systems helps avoid testing the solver against only one specific form of $f(t,x)$.

---

## Numerical methods

### Forward Euler

The Forward Euler method evaluates the derivative at the beginning of each step:

$$
x_{n+1}=x_n+h f(t_n,x_n).
$$

Here $h$ denotes the actual length of the step; the final step may be shorter than the nominal step size.

For the primary test problem $x'=x$, this becomes

$$
x_{n+1}=x_n+h x_n.
$$

Under standard regularity assumptions, the method has first-order global convergence:

$$
E(h)\approx Ch.
$$

### Classical RK4

The classical fourth-order Runge–Kutta method evaluates the right-hand side at four stages:

$$
\begin{aligned}
k_1 &= f(t_n,x_n),\\
k_2 &= f\left(t_n+\frac h2,\ x_n+\frac h2 k_1\right),\\
k_3 &= f\left(t_n+\frac h2,\ x_n+\frac h2 k_2\right),\\
k_4 &= f(t_n+h,\ x_n+h k_3).
\end{aligned}
$$

The solution is then updated using

$$
x_{n+1}=x_n+\frac h6\left(k_1+2k_2+2k_3+k_4\right).
$$

Here, too, $h$ is the actual length of the step. Under standard regularity assumptions and for sufficiently small steps, RK4 has fourth-order global convergence:

$$
E(h)\approx Ch^4.
$$

RK4 uses the same solver interface as Euler, which makes the methods directly comparable on the same ODE systems.

---

## Verification and testing

Automated tests are implemented separately for Euler and RK4 where the tested property is method-specific.

The test suite currently contains **17 GoogleTest test cases** (9 Euler and 8 RK4), all passing in the latest reported run.

The tests cover:

- integration of the primary test problem;
- reaching the requested final time, including cases where floating-point rounding would otherwise add an extra step;
- different initial conditions and right-hand sides;
- time-dependent right-hand sides;
- numerical accuracy against known exact solutions;
- decreasing error with a smaller step size and observed convergence order;
- integration of a two-dimensional harmonic oscillator;
- rejection of invalid integration parameters;
- correct handling of a zero-length integration interval.

### Convergence experiments

For a method with global error

$$
E(h)\approx Ch^p,
$$

reducing the step size by a factor of two gives approximately

$$
\frac{E(h)}{E(h/2)}\approx 2^p.
$$

The experimental order can therefore be estimated as

$$
p\approx\frac{\log(E(h)/E(h/2))}{\log 2}.
$$

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

## Building

The project uses CMake with Eigen and GoogleTest fetched through CMake's `FetchContent` mechanism.

Configure the project:

```bash
cmake -S . -B build
```

Build:

```bash
cmake --build build
```

---

## Running tests

Run the project's ODE tests (excluding tests registered by the Eigen dependency):

```bash
ctest --test-dir build -R '^(EulerSolverTest|RK4SolverTest)\.' --output-on-failure
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
- forward integration with a shortened final step;
- automated verification with GoogleTest.

Current limitations include:

- Only forward integration is supported (`t_end >= t0`, `h > 0`).
- The entire computed trajectory is stored in memory.
- Adaptive step-size control and stiff ODE methods are not implemented.
- State dimension consistency and state-value finiteness are not fully validated.
- Extremely small time steps and extreme floating-point time scales are not comprehensively handled.
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

The current Euler and RK4 implementations form the numerical foundation for the next stages of development.
