using Revise
using LinearAlgebra, SparseArrays, StaticArrays, Plots
using Datasolver
include("../basic_setup.jl")
#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law

dataset = create_dataset(20, x -> bar_E * x, -strain_limit, strain_limit, noise_magnitude = 0.01 * bar_E)
# # solving
resultsL1 = NLP_solver(
	nonlinear_problem,
	dataset;
	use_L1_norm = true,
	random_init_data = false,
	use_data_bounds = true,
)
resultsL2 = NLP_solver(
	nonlinear_problem,
	dataset;
	use_L1_norm = false,
	random_init_data = false,
	use_data_bounds = true,
)


final_L1 = get_final(resultsL1)
final_L2 = get_final(resultsL2)
@show final_L1.E == final_L2.E
@show final_L1.S == final_L2.S
@show final_L1.e - final_L2.e
@show final_L1.s - final_L2.s

@show final_L1.cost - final_L2.cost
plot_results(resultsL1, dataset = dataset)


@time NLP_solver(
	nonlinear_problem,
	dataset;
	use_L1_norm = false,
	random_init_data = false,
	use_data_bounds = true,
)
