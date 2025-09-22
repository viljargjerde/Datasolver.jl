using Revise
using LinearAlgebra, SparseArrays, StaticArrays, Plots
using Datasolver

#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law

# inputs
bar_area = 0.5;     # [m^2] - cross-sectional area of the bar
dims = 2
alpha = 1.0
bar_len = 1.0;      # [m]   - initial length of the bar
bar_E = 1e3;        # [Pa]  - assumed Young_modulus
num_ele = 8;       # [-]   - number of elements
numDataPts = 201;   # [-]   - number of data points, odd number to ensure zero strain is in the dataset
if dims == 1
	bar_distF = [1.8e2]  # [N]   - constant uniform distributed load
else
	bar_distF = [1.8e2, 1.8e1]     # [N]   - constant uniform distributed load
end


# Generate data 
strain_limit = 2 * norm(bar_distF) * bar_len / (bar_E * bar_area);
if norm(bar_distF) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E / 4 * tanh.(10x), -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E * tanh.(10x), -strain_limit, strain_limit)
SE = hcat(dataset.E, dataset.S)
costFunc_ele = (e, s) -> 0.5 * (dataset.C * e^2 + 1 / dataset.C * s^2);

# node vector
num_node = num_ele + 1;
if dims == 1
	node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
	constrained_dofs = get_constrained_dofs([(1, 1), (num_node, 1)], num_ele, dims)
else
	node_vector = [[x, 0.1x] for x in LinRange(0.0, bar_len, num_node)]
	constrained_dofs = get_constrained_dofs([(1, 1), (1, 2), (num_node, 1), (num_node, 2)], num_ele, dims)
end

# # solving
results = directSolverNonLinearBar(
	node_vector,
	constrained_dofs,
	SE,
	costFunc_ele,
	bar_area,
	bar_distF,
	num_ele,
	dataset.C;
	NR_max_iter = 100,
	NR_num_load_step = 100,
	alpha = alpha,
	random_init_data = false,
	DD_max_iter = 10,
	NR_damping = 0.5,
);

plot(results.u)
plot_results(results, dataset = dataset)

