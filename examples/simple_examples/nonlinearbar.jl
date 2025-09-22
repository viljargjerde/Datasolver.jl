using Revise
include("../basic_setup.jl")
using LinearAlgebra, SparseArrays, StaticArrays, Plots
using Datasolver

#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law

# # inputs
# bar_length = 1.0 * 1000      # [mm]   - initial length of the bar
# area = 400.0     # [mm^2] - cross-sectional area of the bar
# force = x -> [0.1]  # [kN/mm]   - constant uniform distributed load
# num_ele = 32       # [-]   - number of elements
# alpha = 1.0

# bar_E = 1e3;        # [MPa]  - Young_modulus
# numDataPts = 128;   # [-]   - number of data points, odd number to ensure zero strain is in the dataset

# problem = fixedBarproblem1D(
# 	bar_length,
# 	area,
# 	force,
# 	num_ele,
# 	alpha;
# 	right_fixed = true,
# )


# # Generate data 
# dist_F = problem.force(nothing)
# strain_limit = norm(dist_F) * problem.length / (2 * bar_E * problem.area);
# if norm(problem.force(1)) <= 1e-6
# 	strain_limit += 1.0
# end
# dataset = create_dataset(numDataPts, x -> bar_E / 4 * tanh.(10x), -strain_limit, strain_limit)
# # dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
# # dataset = create_dataset(numDataPts, x -> bar_E * tanh.(10x), -strain_limit, strain_limit)

# # node vector
# num_node = num_ele + 1;
# node_vector = [[x] for x in LinRange(0.0, problem.length, num_node)]
# constrained_dofs = get_constrained_dofs([(1, 1), (num_node, 1)], problem.num_ele, problem.dims)

# # solving
# results = directSolverNonLinearBar(
# @profview_allocs directSolverNonLinearBar(
# dataset = create_dataset(17, x -> bar_E * x, -strain_limit, strain_limit, noise_magnitude = 0.2)


dataset = create_dataset(2^10, x -> bar_E * x, 0.0, 2 * strain_limit)
# linear_problem, nonlinear_problem = get_problems(6)
nonlinear_problem = fixedBarproblem1D(
	bar_length,
	area,
	x -> [100.0],  # [N/mm]   - constant uniform distributed load
	# x -> [x < 0.9 * L_0 ? 500e3 : 0],  # [N/mm]   - constant uniform distributed load
	# x -> [500e3 / 200 ],  # [N/mm]   - constant uniform distributed load
	512,
	1.0;
	right_fixed = false,
)

@profview Datasolver.greedyLocalSearchSolverNonLinearBar(
	nonlinear_problem,
	dataset;
	random_init_data = false,
	NR_tol = 1e-8,
	cache_ADM = true,
);

# results = directSolverNonLinearBar(
# 	nonlinear_problem,
# 	dataset;
# 	random_init_data = false,
# 	verbose = true,
# 	NR_max_iter = 100,
# 	NR_tol = 1e-9,
# )

@time directSolverNonLinearBar(
	nonlinear_problem,
	dataset;
	random_init_data = false,
	verbose = false,
	NR_max_iter = 100,
	NR_tol = 1e-8,
);
# results.cost

# plot(results.u)
@show norm(results.compatibility)
@show norm(results.equilibrium);
results.cost
plot_results(results, dataset = dataset)
scatter(dataset.E, dataset.S)



