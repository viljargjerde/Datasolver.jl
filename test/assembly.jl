using Test, SparseArrays, LinearAlgebra, Statistics

using Datasolver.DataDrivenNonlinearBar
using Datasolver

@testset "assembly linearized system matrix and rhs" begin
	# inputs
	bar_len = 1.5      # [m]   - initial length of the bar
	bar_area = 2e-3    # [m^2] - cross-sectional area of the bar
	bar_distF = 1.8    # [N]   - constant uniform distributed load
	bar_E = 1e4        # [Pa]  - assumed Young_modulus
	num_ele = 16       # [-]   - number of elements
	numDataPts = 200   # [-]   - number of data points


	# generate data: linear function
	strain_limit = 2 * bar_distF * bar_len / (bar_E * bar_area)
	dataset = create_dataset(numDataPts, x -> bar_E * tanh.(500 .* x), -strain_limit, strain_limit)
	SE = hcat(dataset.E, dataset.S)
	# SE = generateDataHookLaw(Young_modulus = bar_E, numDataPts = numDataPts, strainRange = [-strain_limit, strain_limit])
	# costFunc_constant = mean(SE[:, 2] ./ SE[:, 1])



	# node vector
	num_node = num_ele + 1
	node_vector = collect(LinRange(0.0, bar_len, num_node))


	# data initialization
	init_data_id = rand(1:numDataPts, num_ele)
	init_data = SE[init_data_id, :]
	data_star = deepcopy(init_data)


	# ndofs
	ndof_u = ndof_lambda = num_node
	ndof_e = ndof_s = ndof_mu = num_ele

	ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
	ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]


	# linear system matrix and rhs for the linear case
	A = assembleLinearSystemMatrix(node_vector, num_ele, ndofs, dataset.C, bar_area)

	rhs_lin = assembleRhsLinearBar(data_star = data_star, node_vector = node_vector, num_ele = num_ele, ndofs = ndofs, costFunc_constant = dataset.C, bar_distF = bar_distF, cross_section_area = bar_area)

	# initial guess
	x = zeros(ndof_tot)

	# assembly linearized system and rhs
	rhs = assembleEquilibriumResidual(x = x, data_star = data_star, node_vector = node_vector, num_ele = num_ele, ndofs = ndofs, costFunc_constant = dataset.C, bar_distF = bar_distF, cross_section_area = bar_area)

	J = assembleLinearizedSystemMatrix(x = x, node_vector = node_vector, num_ele = num_ele, ndofs = ndofs, costFunc_constant = dataset.C, cross_section_area = bar_area)

	# tests
	@test norm(J - A) ≈ 0
	@test norm(rhs - rhs_lin) ≈ 0
end
