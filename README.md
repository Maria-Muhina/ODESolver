# Numerical Dynamics in C++

A small C++ project for implementing and testing numerical methods for ordinary differential equations (ODEs).

The project currently focuses on two explicit time-stepping methods: the Forward Euler method and the classical fourth-order Runge–Kutta method (RK4).

The implementation is intentionally small. The goal at this stage is to establish a clean numerical-solver foundation before adding more advanced methods and engineering infrastructure.

---

## Overview

The project considers scalar initial-value problems of the form

\[
x'(t) = f(t, x), \qquad x(t_0) = x_0.
\]

A numerical solver takes an ODE system, an initial condition, a final time, and a nominal step size, and returns the computed solution at the integration points.

The main components are:

- `OdeSystem` — interface for the right-hand side of an ODE;
- `EulerSolver` — Forward Euler implementation;
- `RK4Solver` — classical fourth-order Runge–Kutta implementation;
- concrete ODE systems used for testing;
- GoogleTest-based automated verification.

The solvers use the same basic interface so that different numerical methods can be applied to the same initial-value problem.

---

## Project structure

The main numerical components are organized as follows:

```text
.
├── OdeSystem.h / .cpp
├── EulerSolver.h / .cpp
├── RK4Solver.h / .cpp
├── ExponentialGrowth.h / .cpp
├── DoubleGrowth.h / .cpp
├── TimeDependentSystem.h / .cpp
├── EulerSolverTests.cpp
├── RK4SolverTests.cpp
├── CMakeLists.txt
└── .gitignore
```

`main.cpp` and `Counter.*` are remnants of the initial project skeleton and are not part of the current numerical-solver design.

Generated build and test artifacts are kept outside the source structure and are not part of the project source code.

---

## Mathematical foundation

### Problem formulation

The solver works with an initial-value problem

\[
x'(t) = f(t,x), \qquad x(t_0)=x_0.
\]

Here:

- \(t\) is the independent variable;
- \(x(t)\) is the state;
- \(f(t,x)\) is the right-hand side of the differential equation;
- \(t_0\) and \(x_0\) define the initial condition.

The numerical solution is computed at a sequence of time points

\[
t_0,\ t_1,\ \ldots,\ t_N.
\]

The nominal step size is \(h\). If the requested final time is not reached exactly by a full step, the final step is shortened to end at `t_end`.

### Primary test problem

The main test problem used throughout the project is

\[
x'(t)=x(t), \qquad x(0)=1.
\]

Its exact solution is

\[
x(t)=e^t.
\]

This problem is useful because the exact solution is known, allowing the numerical solution to be compared directly with the analytical one.

It is also used for convergence experiments for both Euler and RK4.

### Additional test systems

The project also contains several other ODE systems:

- \(x'(t)=2x(t)\), used to verify that the solvers work with a different growth rate;
- \(x'(t)=t\), \(x(0)=1\), whose exact solution is
  \[
  x(t)=1+\frac{t^2}{2},
  \]
  used in particular to verify correct handling of time-dependent right-hand sides.

Using several systems helps avoid testing the solver against only one specific form of \(f(t,x)\).

---

## Numerical methods

### Forward Euler

The Forward Euler method evaluates the derivative at the beginning of each step:

\[
x_{n+1}
=
x_n+h f(t_n,x_n).
\]

The implementation also handles a shortened final step when `t_end` is not an exact multiple of the nominal step size.

For the primary test problem \(x'=x\), this becomes

\[
x_{n+1}=x_n+h x_n.
\]

The method has first-order global convergence:

\[
E(h)\approx C h.
\]

### Classical RK4

The classical fourth-order Runge–Kutta method evaluates the right-hand side at four stages:

\[
k_1=f(t_n,x_n),
\]

\[
k_2=f\left(t_n+\frac h2,\,
x_n+\frac h2 k_1\right),
\]

\[
k_3=f\left(t_n+\frac h2,\,
x_n+\frac h2 k_2\right),
\]

\[
k_4=f(t_n+h,\,
x_n+h k_3).
\]

The solution is then updated using

\[
x_{n+1}
=
x_n+
\frac h6
\left(k_1+2k_2+2k_3+k_4\right).
\]

For sufficiently small step sizes, the method has fourth-order global convergence:

\[
E(h)\approx C h^4.
\]

The RK4 implementation uses the same solver interface as Euler, which makes the methods directly comparable on the same ODE systems.

---

## Verification and testing

Automated tests are implemented separately for the Euler and RK4 solvers where the tested property is method-specific.

The test suite currently contains **13 GoogleTest test cases**, all passing.

The tests cover several types of behavior:

- correct integration of the primary test problem;
- correct final time;
- different initial conditions;
- different right-hand sides;
- time-dependent right-hand sides;
- numerical accuracy against known exact solutions;
- reduction of error when the step size is reduced;
- observed convergence order.

### Convergence experiments

For a method with global error

\[
E(h)\approx C h^p,
\]

reducing the step size by a factor of two gives approximately

\[
\frac{E(h)}{E(h/2)}\approx 2^p.
\]

The experimental order can therefore be estimated as

\[
p\approx
\frac{\log(E(h)/E(h/2))}{\log 2}.
\]

#### Euler

For the primary test problem \(x'(t)=x(t)\), \(x(0)=1\), the endpoint error was evaluated at \(t=0.3\) for the following step sizes:

```text
h = 0.1
h = 0.05
h = 0.025
```

The measured absolute errors were approximately:

```text
h       E(h)
0.1     9.76e-3
0.05    4.97e-3
0.025   2.52e-3
```

The corresponding experimental orders were approximately:

```text
0.1  → 0.05   : p ≈ 0.95
0.05 → 0.025  : p ≈ 0.97
```

The estimates approach the expected first-order behavior as the step size decreases.

#### RK4

For the same primary test problem, the endpoint error was evaluated at \(t=0.3\) for:

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
0.2  → 0.1   : p ≈ 3.34
0.1  → 0.05  : p ≈ 3.94
0.05 → 0.025 : p ≈ 3.97
```

The first estimate is lower because the coarser step size is not yet fully in the asymptotic convergence regime. As the step size decreases, the measured order approaches the expected value of 4.

For \(h=0.2\), the interval \([0,0.3]\) is not covered by an integer number of full steps, so the implementation uses a shortened final step. This is part of the solver's normal final-step handling.

The observed increase of these particular order estimates is a property of this experiment; it is not assumed to be a general requirement for every RK4 convergence experiment.

---

## Building

The project uses CMake.

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

Run the complete test suite:

```bash
ctest --test-dir build
```

The test executable can also be run directly:

```bash
./build/ODESolverTests
```

A specific GoogleTest test can be selected with a filter, for example:

```bash
./build/ODESolverTests --gtest_filter=RK4SolverTest.SolvesDoubleGrowth
```

---

## Current limitations

The current implementation is deliberately limited in scope.

At this stage it supports:

- scalar ODEs;
- fixed-step explicit Euler;
- fixed-step classical RK4;
- known scalar test systems;
- automated verification with GoogleTest.

The following are not implemented yet:

- adaptive step-size control;
- stiff ODE methods;
- vector-valued state variables;
- matrix/vector-based numerical linear algebra;
- systematic performance benchmarking;
- profiling;
- continuous integration;
- extended public API documentation.

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