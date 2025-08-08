# 2D Finite Element Poisson Solver with Adaptive Mesh Refinement

**refine2d-fem** is a small, self-contained project that solves the Poisson equation on a 2D domain using the Finite Element Method (FEM) with linear triangular (P1) elements, and includes a starter implementation of Adaptive Mesh Refinement (AMR) using a Zienkiewicz–Zhu (ZZ) error estimator.

It’s written in Python for clarity, but the structure is designed to port cleanly to C++/Eigen for higher performance.

## What This Project Demonstrates
- **Numerical PDE solving** with FEM
- **Adaptive mesh refinement** based on error estimation
- **Clean project architecture**: mesh generation → assembly → solve → error estimate → refinement → export
- **Scientific visualization** with ParaView (`.vtu` output)

## What Are We Solving?

### 1. The Poisson Equation
We solve the 2D Poisson problem:
$
-\nabla \cdot (\kappa \nabla u) = f \quad \text{in } \Omega, \quad u = g \ \text{on } \partial \Omega
$
- $u(x,y)$ = unknown scalar field (think: temperature on a metal plate)
- $f(x,y)$ = source term (heat generation)
- $kappa(x,y)$ = material conductivity
- $g(x,y)$ = fixed boundary values

**Example physical analogy:**  
Imagine a thin square metal plate. You heat certain spots (via $f$) and keep the edges fixed at certain temperatures (via $g$). We want to know the temperature everywhere in the plate.

### 2.  The Finite Element Method (FEM)
1.  Break the domain into small shapes (triangles).
2. Assume a simple function (linear in our case) on each triangle.
3. Stitch them together so they agree at the edges.
4. Write equations saying the approximate solution should satisfy the PDE in an average sense (weak form).
5. Solve the resulting system of linear equations for the values at the triangle corners (nodes).

### 3. Adaptive Mesh Refinement
Instead of refining everywhere, we:
1. Estimate the error in each element (ZZ error estimator).
2. Mark elements where the error is largest.
3. Refine only those elements, creating a denser mesh where needed.

**Why?** This saves computation and you get accuracy where the solution changes rapidly, without wasting work where the solution is smooth.

## The Math in Short

### 1. Weak form
Multiply PDE by a test function $v$, integrate, integrate by parts:
$$
a(u,v) = \int_{\Omega} \kappa \nabla u \cdot \nabla v \, dx, \quad
l(v) = \int_{\Omega} f v \, dx
$$
Find $u_h$ in the finite element space $V_h$ such that:
$$
a(u_h, v_h) = l(v_h) \quad \forall v_h \in V_h
$$

### 2. Linear triangular (P1) basis
Each basis function is 1 at its node, 0 at others, linear on each triangle.

### 3. Element stiffness
$$
K^e_{ij} = \int_{T_e} \kappa \, \nabla \phi_i \cdot \nabla \phi_j \, dA
$$

### 4. Load vector
$$
f^e_i = \int_{T_e} f \, \phi_i \, dA
$$

### 5. ZZ error estimator
- Compute gradient of $u_h$ on each triangle.
- Recover a smoother “nodal gradient” by averaging over neighbors.
- Compare the two; large differences mean large error.

## Quick Start
1. Clone or download
```bash
git clone https://github.com/harmanxgill/refine2d-fem.git
cd refine2d-fem
```

2. Install dependencies
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

3. Run the example
```bash
python examples/run_poisson.py
```

This writes `out/solution_cycle0.vtu`; open it in ParaView to see the solution and error indicator.

##  Next Steps & Roadmap
- Implement newest-vertex bisection for conforming refinement.
- Implement Dörfler marking (bulk criterion) instead of top-fraction marking.
- Manufactured solution tests with H1 and L2 error convergence plots.
- Optional: swap SciPy CG solver for PyAMG preconditioned CG.
- C++/Eigen port for performance.

## References
- Zienkiewicz, O. C., Taylor, R. L., & Zhu, J. Z. (2013). The finite element method: Its basis and fundamentals (7th ed.). Butterworth-Heinemann. https://doi.org/10.1016/C2009-0-24909-9
- Hughes, T. J. R. (2000). The finite element method: Linear static and dynamic finite element analysis. Dover Publications. (Original work published 1987)
- Ladevèze, P., & Pelle, J.-P. (2004). Mastering calculation in linear and nonlinear mechanics. Springer. https://doi.org/10.1007/978-0-387-21787-2
- Ainsworth, M., & Oden, J. T. (1997). A posteriori error estimation in finite element analysis. Computer Methods in Applied Mechanics and Engineering, 142(1–2), 1–88. https://doi.org/10.1016/S0045-7825(96)01107-3
