# Mixed-Precision Iterative Refinement with Stochastic Rounding

This repository contains the code, numerical experiments, and selected plots for my honors thesis on **mixed-precision iterative refinement (IR)** and **stochastic rounding (SR)** in ill-conditioned problems.

The project studies whether stochastic rounding in low-precision computations can help mitigate stagnation and maintain accuracy in mixed-precision iterative refinement algorithms. The experiments focus on inverse problems and image deblurring examples, including the **inverse heat problem** and the **PRblur image restoration problem**.

## Project Overview

Mixed-precision algorithms use lower precision arithmetic to reduce computational cost while relying on higher precision steps to preserve accuracy. In this project, stochastic rounding is used in the low-precision computation step of iterative refinement and related solvers.

Main topics explored in this repository include:
- iterative refinement for linear systems,
- LU-based and GMRES-based inner solves,
- stochastic rounding versus round-to-nearest,
- Tikhonov regularization for ill-posed problems,
- behavior under different noise levels and precisions.

## Repository Contents

Typical contents of this repository include:

- `src/` — core Julia implementations of the algorithms  
- `experiments/` — scripts for running numerical tests  
- `plots/` — saved figures from experiments  
- `data/` — problem data such as PRblur test data  
- `results/` — output files, relative error histories, or summaries  
- `README.md` — project description and usage instructions  

You may rename these folders to match the actual repository structure.

## Numerical Problems

The main test problems in this repository are:

### 1. Inverse Heat Problem
A severely ill-conditioned inverse problem used to test the behavior of iterative refinement and regularization.

### 2. PRblur Image Deblurring Problem
A large-scale image restoration problem with Kronecker product structure, used to compare regularized and unregularized behavior under different precisions and noise levels.

## Methods Implemented

This repository includes code related to:
- standard iterative refinement,
- mixed-precision iterative refinement,
- LU-based correction solves,
- GMRES-based correction solves,
- stochastic rounding in low precision,
- Tikhonov regularization.

In the experiments, the residual is typically computed in higher precision, while the correction step is solved in lower precision.

## Software and Packages

The code is written in **Julia**.

Packages used in this project may include:
- `LinearAlgebra`
- `Plots`
- `IterativeSolvers`
- `StochasticRounding`
- `MAT`

Install dependencies in Julia with:

```julia
using Pkg
Pkg.add(["LinearAlgebra", "Plots", "IterativeSolvers", "StochasticRounding", "MAT"])
