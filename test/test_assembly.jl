using Test, SparseArrays, LinearAlgebra, Statistics


@testset "assembly linearized system matrix and rhs" begin
	# inputs
	bar_length = 1.5      # [m]   - initial length of the bar
	area = 2e-3    # [m^2] - cross-sectional area of the bar
	force = x -> [1.8]    # [N]   - constant uniform distributed load
	num_ele = 16       # [-]   - number of elements
	numDataPts = 200   # [-]   - number of data points
	problem_lin = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		0.0;
		right_fixed = true,
	)

	problem_nonlin = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		1.0;
		right_fixed = true,
	)

	bar_E = 1e4        # [Pa]  - assumed Young_modulus
	# generate data: linear function
	strain_limit = 2 * norm(force(0)) * bar_length / (bar_E * area)
	dataset = create_dataset(numDataPts, x -> bar_E * tanh.(500 .* x), -strain_limit, strain_limit)

	num_node = num_ele + 1
	# data initialization
	init_data_id = rand(1:numDataPts, num_ele)
	E, S = dataset[init_data_id]

	# initial guess
	x = zeros(sum(Datasolver.get_ndofs(problem_nonlin)))
	# linear system matrix and rhs for the linear case

	A = Datasolver.assembleLinearizedSystemMatrix(x, problem_lin, dataset.C)

	rhs_lin = Datasolver.assembleEquilibriumResidual(x, E, S, dataset.C, problem_lin, 1.0)



	J = Datasolver.assembleLinearizedSystemMatrix(x, problem_nonlin, dataset.C)

	rhs = Datasolver.assembleEquilibriumResidual(x, E, S, dataset.C, problem_nonlin, 1.0)

	# tests
	@test norm(J - A) ≈ 0
	@test norm(rhs - rhs_lin) ≈ 0
end
