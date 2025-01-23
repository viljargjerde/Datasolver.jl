using LinearAlgebra, SparseArrays, StaticArrays, Statistics
using Datasolver


#### data-driven for linear fixed-free bar with constant uniform distributed load and linear material law


# inputs
bar_len = 2.0;      # [m]   - initial length of the bar
bar_area = 0.5;     # [m^2] - cross-sectional area of the bar
bar_distF = 1.8e2;    # [N]   - constant uniform distributed load
bar_E = 1e4;        # [Pa]  - assumed Young_modulus
num_ele = 10;       # [-]   - number of elements
numDataPts = 200;   # [-]   - number of data points


strain_limit = 2 * bar_distF * bar_len / (bar_E * bar_area);
dataset = dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
SE = hcat(dataset.E, dataset.S)

# elementwise cost function
costFunc_ele = (e, s) -> 0.5 * (dataset.C * e^2 + 1 / dataset.C * s^2);


# node vector
num_node = num_ele + 1;
node_vector = collect(LinRange(0.0, bar_len, num_node));

# solving
results = Datasolver.DataDrivenNonlinearBar.directSolverLinearBar(node_vector = node_vector, data_set = SE, costFunc_ele = costFunc_ele, num_ele = num_ele, costFunc_constant = dataset.C, bar_distF = bar_distF, cross_section_area = bar_area);

# plots 

Datasolver.plot_results(results, dataset = dataset)
