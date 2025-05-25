# DataSolver: A Data-Driven Solver Package

**Author**: Viljar Helgestad Gjerde  
**Project**: A Data-Driven Model for the Analysis of Geometrically Nonlinear One-Dimensional Structures

## Overview

`Datasolver.jl` is a Julia package developed as part of a master's thesis, focusing on solving data-driven boundary value problems using Data-Driven Computational Mechanics ([DDCM](https://www.sciencedirect.com/science/article/abs/pii/S0045782516300238)). Instead of relying on traditional constitutive models, this DDCM matches local strain-stress states directly against experimental or synthetic datasets. The package includes three solvers, all supporting both linear and nonlinear strain measures. For more information, refer to the master thesis "A Data-Driven Model for the Analysis of Geometrically Nonlinear One-Dimensional Structures" that is to be published by the University of Bergen.

This package is designed for research purposes, facilitating experimentation, visualization, and benchmarking in academic settings.

---

## Key Features

-  **Data-Driven Formulation**: Solve mechanical problems using empirical data rather than constitutive laws.
-  **Modular Problem Setup**: Define 1D bar problems with various loading conditions and boundary constraints.
-  **Multiple Solvers**: Includes a mixed-integer NLP solver and a greedy search solver for bar elements.
-  **Postprocessing & Visualization**: Generate convergence plots, stress-strain diagrams, and more.
-  **Convergence Analysis**: Benchmark accuracy across datasets and discretization levels.
-  **Custom Dataset Tools**: Create synthetic datasets with optional noise.

---

## Installation

This package is not registered in the Julia package registry. To install locally:

```julia
import Pkg
Pkg.add(url = "https://github.com/viljargjerde/Datasolver.jl")
```

Ensure the following dependencies are installed:

```julia
Pkg.add(["Gurobi", "JuMP", "Plots", "JSON", "DataFrames", "Dierckx", "StatsBase"])
```

Note: Proper configuration of [Gurobi](https://www.gurobi.com/) is required for the MINLP-based solver.

---

## Example Usage

### 1. Generate a Synthetic Dataset

```julia
using Datasolver

f = x -> 5e6 * tanh(500 * x)
dataset = create_dataset(100, f, -0.005, 0.005, noise_magnitude = 0.01)
plot_dataset(dataset)
```

### 2. Define a 1D Bar Problem

```julia
bar = fixedBarproblem1D(length = 1.0, area = 1.0, force = x -> [100.0], num_ele = 10, alpha = 0.0, right_fixed = true)
```

### 3. Solve with NLP or Greedy Solver

```julia
result = NLP_solver(bar, dataset, use_L1_norm = true) # MINLP
# or
result = directSolverNonLinearBar(bar, dataset) # ADM
# or
result = greedyLocalSearchSolverNonLinearBar(bar, dataset) # GO-ADM
```

### 4. Plot Results

```julia
plot_results(result, dataset = dataset, title = "Solver Output")
```

---




## Project Structure

- `src/` – Core module files
  - `Datasolver.jl` – Main module entry point
  - `dataset.jl` – Dataset structures and utilities
  - `barproblem.jl` – Bar problem formulations
  - `LP_solver.jl` – MINLP solver using JuMP and Gurobi
  - `utils.jl` – Visualization and analysis tools
  - `assembly.jl` – Core numerical assembly routines
- `examples/` – Scripts for all examples presented in the thesis
- `test/` – Unit tests for the package

---
