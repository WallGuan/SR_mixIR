{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Mixed-Precision Iterative Refinement with Stochastic Rounding\
\
This repository contains the code, numerical experiments, and selected plots for my honors thesis on **mixed-precision iterative refinement (IR)** and **stochastic rounding (SR)** in ill-conditioned problems.\
\
The project studies whether stochastic rounding in low-precision computations can help mitigate stagnation and maintain accuracy in mixed-precision iterative refinement algorithms. The experiments focus on inverse problems and image deblurring examples, including the **inverse heat problem** and the **PRblur image restoration problem**.\
\
## Project Overview\
\
Mixed-precision algorithms use lower precision arithmetic to reduce computational cost while relying on higher precision steps to preserve accuracy. In this project, stochastic rounding is used in the low-precision computation step of iterative refinement and related solvers.\
\
Main topics explored in this repository include:\
- iterative refinement for linear systems,\
- LU-based and GMRES-based inner solves,\
- stochastic rounding versus round-to-nearest,\
- Tikhonov regularization for ill-posed problems,\
- behavior under different noise levels and precisions.\
\
## Repository Contents\
\
Typical contents of this repository include:\
\
- `src/` \'97 core Julia implementations of the algorithms  \
- `experiments/` \'97 scripts for running numerical tests  \
- `plots/` \'97 saved figures from experiments  \
- `data/` \'97 problem data such as PRblur test data  \
- `results/` \'97 output files, relative error histories, or summaries  \
- `README.md` \'97 project description and usage instructions  \
\
You may rename these folders to match the actual repository structure.\
\
## Numerical Problems\
\
The main test problems in this repository are:\
\
### 1. Inverse Heat Problem\
A severely ill-conditioned inverse problem used to test the behavior of iterative refinement and regularization.\
\
### 2. PRblur Image Deblurring Problem\
A large-scale image restoration problem with Kronecker product structure, used to compare regularized and unregularized behavior under different precisions and noise levels.\
\
## Methods Implemented\
\
This repository includes code related to:\
- standard iterative refinement,\
- mixed-precision iterative refinement,\
- LU-based correction solves,\
- GMRES-based correction solves,\
- stochastic rounding in low precision,\
- Tikhonov regularization.\
\
In the experiments, the residual is typically computed in higher precision, while the correction step is solved in lower precision.\
\
## Software and Packages\
\
The code is written in **Julia**.\
\
Packages used in this project may include:\
- `LinearAlgebra`\
- `Plots`\
- `IterativeSolvers`\
- `StochasticRounding`\
- `MAT`\
\
Install dependencies in Julia with:\
\
```julia\
using Pkg\
Pkg.add(["LinearAlgebra", "Plots", "IterativeSolvers", "StochasticRounding", "MAT"])}