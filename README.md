# DataSolver: A Data-Driven Solver Package

**Author: Viljar Helgestad Gjerde**  
**Date: 04.10.24**

## Overview

DataSolver is a Julia package that provides an easy-to-use interface for setting up and solving data-driven engineering problems. It is developed as part of Viljar Helgestad Gjerde's master's thesis. The package includes tools for constructing synthetic data driven problems and solving them with the included solver. The code is under active development and will change with no warning. 

### Features
- **Data-driven Solver** (`datasolve`): Solves a data driven problem with extensive logging for each iteration. 
- **Flexible Dataset Creation**: Easily generate synthetic datasets for strain-stress analysis.
- **Problem Setup Functions**: Quickly define nodes, elements, and connections for structural problems.
- **Visualization Tools**: Visualize initial configurations, datasets, solver results, and perform convergence analysis.
- **Example Notebook**: Contains an interactive Pluto notebook for demonstrating the use of the package. 

## Prerequisites

Make sure you have Julia installed on your machine. You can get Julia from [https://julialang.org/downloads/](https://julialang.org/downloads/).

To run the Pluto.jl notebook that demonstrates the package, you will need to install `Pluto.jl` as well as `PlutoUI`. You can install these by running:

```julia
using Pkg
Pkg.add(["Pluto", "PlutoUI"])
```

## Getting Started

### Loading the Package

First, clone the repository. You can then add the package by running 
```julia
] add ./Datasolver
```
or
```julia
import Pkg; Pkg.add(path = "./Datasolver")
```

If you are working in another directory, the paths will need to change to reflect this.

### Example Notebook

The package includes a Pluto notebook (`DataSolverDemo.jl`) that demonstrates the usage of the DataSolver package. You can open it with:

```julia
using Pluto
Pluto.run()
```

Navigate in the browser to the location where the notebook/repository is saved and open it to explore the interactive demonstration. Note that some of the cells in the notebook are hidden for presentation purposes, but can be shown by clicking the eye icon next to the cell. 

## Usage Guide

### Creating a Dataset

The `create_dataset` function allows you to generate a synthetic dataset for strain-stress analysis:

```julia
N = 100
min_strain = -0.005
max_strain = 0.005
strain_stress_relation = x -> 5e6 * tanh(500 * x)
noise_magnitude = 1e4

# Create dataset
my_dataset = create_dataset(N, strain_stress_relation, min_strain, max_strain, noise_magnitude=noise_magnitude)
```
This creates a dataset containing `N` data points, simulating a strain-stress relation defined by `strain_stress_relation`.

### Setting up and Solving a Problem

You can set up a structural problem using the `setup_1d_bar_problem` function. This function defines the nodes, elements, and discretized forces for a bar problem:

```julia
N_elements = 20
L = 3.6
force_function = x -> 4000 * sin(pi * x / L)

connections, Φ, f, fixed_dofs = setup_1d_bar_problem(N_elements, L, force_function)
```

The `datasolve` or `my_new_solver` function can then be used to solve the problem:

```julia
result = datasolve(connections, Φ, 0.002, my_dataset, f, fixed_dofs, verbose=false)
```

### Visualizing the Results

DataSolver provides several functions for visualizing the results, including `plot_dataset`, `plot_configuration`, and `plot_results`:

```julia
# Plot the initial configuration of nodes and connections
plot_configuration(Φ, connections)

# Plot the dataset
plot_dataset(my_dataset)

# Plot the dataset with the chosen points from the result
plot_dataset(my_dataset, get_final(result))

# Plot the results of the solver
plot_results(result, dataset=my_dataset)
```

### Convergence Analysis

DataSolver allows you to perform convergence analysis to see how different numbers of elements and data points affect the accuracy of the solution:

```julia
convergence_analysis(results, analytical_u)
```
This function generates a contour plot representing the relative differences between the solved and analytical displacements.

## Notebook Demo

The included Pluto notebook demonstrates several examples of using the `DataSolver` package, including:
- Setting up different constitutive relationships.
- Performing convergence analysis.
- Comparing different solver methods.
- Visualization of strain-stress datasets and solver results.

You can adjust parameters such as noise, number of data points, and elements interactively using the sliders provided in the Pluto notebook.


[![Build Status](https://github.com/viljargjerde/Datasolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/viljargjerde/Datasolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)
