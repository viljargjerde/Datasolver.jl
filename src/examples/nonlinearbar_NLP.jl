using Revise
using LinearAlgebra, SparseArrays, StaticArrays, Statistics, Plots
using Datasolver


#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law

# inputs
bar_area = 0.5;     # [m^2] - cross-sectional area of the bar
alpha = 1.0
bar_len = 1.0;      # [m]   - initial length of the bar
bar_E = 1e3;        # [Pa]  - assumed Young_modulus
num_ele = 4;       # [-]   - number of elements
numDataPts = 20;   # [-]   - number of data points, odd number to ensure zero strain is in the dataset
bar_distF = 1.8e4  # [N]   - constant uniform distributed load

# generate data: linear function 

strain_limit = 0.5 * norm(bar_distF) * bar_len / (bar_E * bar_area);
if norm(bar_distF) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E / 4 * tanh.(10x), -strain_limit, strain_limit)
@show dataset.C
SE = hcat(dataset.E, dataset.S)

# node vector
num_node = num_ele + 1;
# node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
node_vector = [x for x in LinRange(0.0, bar_len, num_node)]
constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1), (num_node, 1)], num_ele, 1)

# # solving
results = Datasolver.nonlin_LP_solver(
	node_vector,
	dataset,
	bar_area,
	bar_distF,
	constrained_dofs,
)

plot(results.u, marker = :x, label = "u")
Datasolver.plot_results(results, dataset = dataset)
# plot(1:length(costFunc_global), costFunc_global, xscale = :log10, yscale = :log10)

# results.u
