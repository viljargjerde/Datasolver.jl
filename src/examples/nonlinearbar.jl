using Revise
using LinearAlgebra, SparseArrays, StaticArrays, Statistics, Plots
using Datasolver


#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law


# inputs
bar_area = 0.5;     # [m^2] - cross-sectional area of the bar
dims = 1
alpha = 1.0

bar_len = 9.0;      # [m]   - initial length of the bar
bar_E = 1e4;        # [Pa]  - assumed Young_modulus
num_ele = 10;       # [-]   - number of elements
numDataPts = 200;   # [-]   - number of data points
if dims == 1
	bar_distF = [1.8e2]    # [N]   - constant uniform distributed load
else
	bar_distF = [1.8e2, 1.8e1]    # [N]   - constant uniform distributed load
end
# bar_len = 2.0;      # [m]   - initial length of the bar
# bar_area = 1.5;     # [m^2] - cross-sectional area of the bar
# bar_distF = 1.8;    # [N]   - constant uniform distributed load
# bar_E = 1e4;        # [Pa]  - assumed Young_modulus
# num_ele = 10;       # [-]   - number of elements
# numDataPts = 200;   # [-]   - number of data points

# generate data: linear function
strain_limit = 2 * norm(bar_distF) * bar_len / (bar_E * bar_area);
dataset = dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
SE = hcat(dataset.E, dataset.S)
costFunc_ele = (e, s) -> 0.5 * (dataset.C * e^2 + 1 / dataset.C * s^2);

# node vector
num_node = num_ele + 1;
# node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
if dims == 1
	node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
	constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1)], num_ele, dims)
else
	node_vector = [[x, 0.1x] for x in LinRange(0.0, bar_len, num_node)]
	constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1), (1, 2)], num_ele, dims)
end
# constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1)], num_ele, length(node_vector[1]))

# # solving
results = Datasolver.DataDrivenNonlinearBar.directSolverNonLinearBar(
	node_vector,
	constrained_dofs,
	SE,
	costFunc_ele,
	bar_area,
	bar_distF,
	num_ele,
	dataset.C;
	NR_max_iter = 100,
	NR_num_load_step = 10,
	alpha = alpha,
);


Datasolver.plot_results(results, dataset = dataset)
# plot(1:length(costFunc_global), costFunc_global, xscale = :log10, yscale = :log10)

