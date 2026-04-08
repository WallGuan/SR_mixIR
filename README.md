# Stochastic Rounding in Mixed-Precision Iterative Refinement
![Julia](https://img.shields.io/badge/Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white)
![MATLAB](https://img.shields.io/badge/MATLAB-E16737?style=for-the-badge)
![Kronecker Product](https://img.shields.io/badge/Kronecker%20Product-0A66C2?style=for-the-badge)
![Iterative Refinement](https://img.shields.io/badge/Iterative%20Refinement-1F6FEB?style=for-the-badge)
![Stochastic Rounding](https://img.shields.io/badge/Stochastic%20Rounding-7A3E9D?style=for-the-badge)
![Mixed Precision](https://img.shields.io/badge/Mixed%20Precision-0E7490?style=for-the-badge)
![GMRES](https://img.shields.io/badge/GMRES-15803D?style=for-the-badge)
![LU](https://img.shields.io/badge/LU-F59E0B?style=for-the-badge)

This repository serves as the GitHub page for the numerical experiments developed for my honors thesis, **Stochastic Rounding in Mixed-Precision Iterative Refinement**. It contains the main Julia codes used in the experiments, supporting Matlab files, output logs and graphs, and related background readings used for the literature review.

The repository is organized mainly around the implementation and testing of mixed-precision iterative refinement algorithms, with stochastic rounding used in the low-precision computation step. The numerical experiments in this project focus primarily on the inverse heat problem and the PRblur image restoration problem.

## Repository Structure

At the top level, the repository contains the following folders:

- **`Main/`**  
  This is the main working directory of the project. Most of the numerical experiment codes are located here.

- **`KroneckerProductTools_matlab/`**  
  This folder contains the original Matlab codes for the Kronecker product tools and object definitions. These files serve as a reference for the Julia implementation.

- **`Other_Julia_code/`**  
  This folder contains additional Julia scripts that are mostly experimental, preliminary, or unrelated to the main numerical experiments reported in the thesis.

- **`Related_papers/`**  
  This folder contains papers and background readings related to the literature review and theoretical foundation of the thesis.

## Main Folder Contents

Most of the codes used for the thesis experiments are inside the `Main/` folder.

### Subfolders

- **`kronecker/`**  
  This folder contains the Julia implementation of the Kronecker product object definition. It is translated and adapted from the Matlab Kronecker product tools.

- **`logs_graphs/`**  
  This folder contains output logs, recorded numerical values, and generated graphs, especially plots of relative error and convergence behavior.

### Dataset

- **`PRblurData_new.mat`**  
  This file contains the data for the PRblur problem, including the Kronecker factors, the right-hand side vector `b`, and the true solution `x_true`.
  
### Main Julia Files

- **`refinement_funcs.jl`**  
  Contains the iterative refinement algorithms and related functions. It also includes the predefined GMRES function used in the experiments.

- **`other_funcs.jl`**  
  Contains supporting helper functions used throughout the experiments to improve efficiency and organization.

- **`heat_problem.jl`**  
  Runs the numerical experiments for the inverse heat matrix problem.

- **`kron_problem.jl`**  
  Runs the numerical experiments for the PRblur problem, which uses Kronecker product structure.

- **`blur_mat.jl`**  
  Outputs the restored solution into `.mat` files so that the image results can be displayed in Matlab.

## Project Overview

This project studies the behavior of **mixed-precision iterative refinement (IR)** with **stochastic rounding (SR)**. The main goal is to examine whether stochastic rounding in low-precision arithmetic can help reduce stagnation and maintain stable convergence in ill-conditioned problems.

The experiments compare different precisions and rounding modes, including standard round-to-nearest and stochastic rounding, under both regularized and unregularized settings.

The main problems considered are:

- **Inverse heat problem**  
  A severely ill-conditioned inverse problem used to test convergence and regularization behavior.

- **PRblur image restoration problem**  
  A large-scale deblurring problem with Kronecker product structure, tested in both noise-free and noisy settings to examine the implementation of SR in a more realistic problem.

## Methods Included

This repository includes implementations related to:

- mixed-precision iterative refinement,
- LU-based correction solves,
- GMRES-based correction solves,
- stochastic rounding in low precision,
- Tikhonov regularization,
- relative error tracking and convergence plots.

In the experiments, the residual and update are typically computed in higher precision, while the correction step is solved in lower precision.

## Outputs

The repository includes numerical outputs such as:

- convergence histories,
- relative error plots, 
- experiment logs,
- restored image outputs saved as `.mat` files. 

In many plots, the first plotted point is **iteration 0**, which represents the relative error of the initial guess. Therefore, a graph may show one additional iteration beyond the number of iterations actually performed by the algorithm.

## Software

The main codes are written in **Julia**, with some supporting and original reference codes in **Matlab**. The code was run on the CPU of an NVIDIA H100 system through a connection to Emory’s ADA server.

The implementation of stochastic rounding in this project uses the Julia package **StochasticRounding.jl**, available at: [StochasticRounding.jl](https://github.com/milankl/StochasticRounding.jl).

## Thesis Context

This repository was developed as part of my honors thesis on stochastic rounding in mixed-precision iterative refinement. It is intended to document the computational experiments, supporting implementations, and outputs used in the thesis work.
