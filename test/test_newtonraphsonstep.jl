using Test, LinearAlgebra, SparseArrays, StaticArrays, Statistics


@testset "NewtonRaphsonStep" begin

	# inputs
	bar_length = 1.5      # [m]   - initial length of the bar
	area = 2e-3    # [m^2] - cross-sectional area of the bar
	force = x -> [1.8]    # [N]   - constant uniform distributed load
	bar_E = 1e4        # [Pa]  - assumed Young_modulus
	num_ele = 10       # [-]   - number of elements
	numDataPts = 200   # [-]   - number of data points

	lin_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		0.0;
		right_fixed = true,
	)

	nonlin_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		1.0;
		right_fixed = true,
	)

	# generate data: linear function
	strain_limit = 2 * norm(force(0)) * bar_length / (bar_E * area)
	dataset = create_dataset(numDataPts, x -> bar_E * tanh.(500 .* x), -strain_limit, strain_limit)

	# randomly assign inital data
	init_data_id = rand(1:numDataPts, num_ele)
	E, S = dataset[init_data_id]


	ndof_tot = sum(Datasolver.get_ndofs(nonlin_problem))


	## solving system with linear strain measure
	# assembly system matrix
	A = Datasolver.assembleLinearizedSystemMatrix(zeros(ndof_tot), lin_problem, dataset.C)

	rhs = Datasolver.assembleEquilibriumResidual(
		zeros(ndof_tot),
		E,
		S,
		dataset.C,
		lin_problem,
	)


	# boundary conditions: fixed-free
	ids = collect(1:size(A, 1))
	deleteat!(ids, lin_problem.constrained_dofs)
	A = A[ids, ids]
	rhs = rhs[ids]


	# solving
	x_lin = qr(Matrix(A)) \ rhs


	## solving system with nonlinear strain measure- test whether the first newton-raphson step leads to the same results
	x_nonlin = Datasolver.NewtonRaphsonStep(
		zeros(ndof_tot),
		E,
		S,
		dataset.C,
		nonlin_problem,
		ids,
		false,
	)

	@test isapprox(norm(x_nonlin[ids] - x_lin), 0, atol = 1e-4, rtol = 1e-4)
end
