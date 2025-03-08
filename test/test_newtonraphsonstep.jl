using Test, LinearAlgebra, SparseArrays, StaticArrays, Statistics


@testset "NewtonRaphsonStep" begin

	# inputs
	bar_len = 1.5      # [m]   - initial length of the bar
	bar_area = 2e-3    # [m^2] - cross-sectional area of the bar
	bar_distF = [1.8]    # [N]   - constant uniform distributed load
	bar_E = 1e4        # [Pa]  - assumed Young_modulus
	num_ele = 10       # [-]   - number of elements
	numDataPts = 200   # [-]   - number of data points


	# generate data: linear function
	strain_limit = 2 * norm(bar_distF) * bar_len / (bar_E * bar_area)
	dataset = create_dataset(numDataPts, x -> bar_E * tanh.(500 .* x), -strain_limit, strain_limit)
	SE = hcat(dataset.E, dataset.S)


	# randomly assign inital data
	init_data_id = rand(1:numDataPts, num_ele)
	init_data = SE[init_data_id, :]
	data_star = deepcopy(init_data)


	# node vector
	num_node = num_ele + 1
	node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]



	# ndofs
	num_node = length(node_vector)
	ndof_u = ndof_lambda = num_node
	ndof_e = ndof_s = ndof_mu = num_ele

	ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
	ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]


	## solving system with linear strain measure
	# assembly system matrix
	A = Datasolver.assembleLinearizedSystemMatrix(zeros(ndof_tot), node_vector, num_ele, ndofs, dataset.C, bar_area, 0.0)

	rhs = Datasolver.assembleEquilibriumResidual(zeros(ndof_tot), data_star, node_vector, num_ele, ndofs, dataset.C, bar_distF, bar_area, 0.0)


	# boundary conditions: fixed-free
	constrained_dofs = [1 (ndof_u + ndof_e + ndof_s + ndof_mu + 1)]
	ids = collect(1:size(A, 1))
	deleteat!(ids, constrained_dofs)
	A = A[ids, ids]
	rhs = rhs[ids]


	# solving
	x_lin = qr(Matrix(A)) \ rhs


	## solving system with nonlinear strain measure- test whether the first newton-raphson step leads to the same results
	x_nonlin = Datasolver.NewtonRaphsonStep(
		zeros(ndof_tot),
		data_star,
		node_vector,
		num_ele,
		ndofs,
		dataset.C,
		bar_distF,
		bar_area,
		1.0,
		constrained_dofs,
	)

	@test norm(x_nonlin - x_lin) ≈ 0
end
