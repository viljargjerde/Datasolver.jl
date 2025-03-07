using Revise
using LinearAlgebra, SparseArrays, StaticArrays, Statistics, Plots
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
use_NR = true
if dims == 1
	bar_distF = [1.8e2]  # [N]   - constant uniform distributed load
else
	bar_distF = [1.8e2, 1.8e1]     # [N]   - constant uniform distributed load
end
# bar_len = 2.0;      # [m]   - initial length of the bar
# bar_area = 1.5;     # [m^2] - cross-sectional area of the bar
# bar_distF = 1.8;    # [N]   - constant uniform distributed load
# bar_E = 1e4;        # [Pa]  - assumed Young_modulus
# num_ele = 10;       # [-]   - number of elements
# numDataPts = 200;   # [-]   - number of data points

# generate data: linear function 

strain_limit = 2 * norm(bar_distF) * bar_len / (bar_E * bar_area);
if norm(bar_distF) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E / 4 * tanh.(10x), -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E * tanh.(10x), -strain_limit, strain_limit)
@show dataset.C
SE = hcat(dataset.E, dataset.S)
costFunc_ele = (e, s) -> 0.5 * (dataset.C * e^2 + 1 / dataset.C * s^2);

# node vector
num_node = num_ele + 1;
# node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
if dims == 1
	node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
	constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1), (num_node, 1)], num_ele, dims)
else
	node_vector = [[x, 0.1x] for x in LinRange(0.0, bar_len, num_node)]
	constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1), (1, 2), (num_node, 1), (num_node, 2)], num_ele, dims)
	# constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1), (1, 2)], num_ele, dims)
end
# constrained_dofs = Datasolver.DataDrivenNonlinearBar.get_constrained_dofs([(1, 1)], num_ele, length(node_vector[1]))



# ####################################################




# ## initialize e_star and s_star
# numDataPts = size(SE, 1)
# num_node = length(node_vector)
# ndof_u = ndof_lambda = num_node * dims
# ndof_e = ndof_s = ndof_mu = num_ele

# ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
# ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]


# init_data = zeros(num_ele, 2)
# s = Datasolver.DataDrivenNonlinearBar.get_initialization_s(bar_distF, node_vector, bar_area, num_ele, ndofs, constrained_dofs, 2)
# S = zeros(length(s))
# E = zeros(length(s))
# for i in eachindex(S)
# 	# choose E, S pair where S is closest to s
# 	best_idx = argmin((abs(SE[j, 2] - s[i]) for j in eachindex(SE[:, 1])))
# 	init_data[i, :] = SE[best_idx, :]
# end


# g_fun = (x, p) -> Datasolver.DataDrivenNonlinearBar.assembleSystemg(x, init_data, node_vector, num_ele, ndofs, dataset.C, bar_area, bar_distF, alpha, constrained_dofs)
# x_init = rand(ndof_tot)
# problem = NonlinearSolve.NonlinearProblem(g_fun, x_init, zeros(ndof_tot))
# sol = solve(problem, NewtonRaphson(), maxiters = 2000)
# ####################################################

# sol.u



###################################
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
	NR_num_load_step = 100,
	alpha = alpha,
	random_init_data = false,
	DD_max_iter = 10,
	NR_damping = 0.5,
	use_NR = use_NR,
);

plot(results.u)
get_final(results)
Datasolver.plot_results(results, dataset = dataset)
# plot(1:length(costFunc_global), costFunc_global, xscale = :log10, yscale = :log10)

