using Test, SparseArrays, LinearAlgebra

using Datasolver.DataDrivenNonlinearBar



@testset "LagrangePolynomials" begin
	# linear Lagrange polynomials
	L0, L1 = linearLagrangePolynomials([0.0, 1.0])

	@test L0(-1) == 1
	@test L0(0) == 0.5
	@test L0(1) == 0

	@test L1(-1) == 0
	@test L1(0) == 0.5
	@test L1(1) == 1


	# 1st derivative of linear Lagrange polynomials
	dL0, dL1 = compute1stDeriv4linearLagrangePolynomials([0.0, 1.0])
	xx = LinRange(-1, 1, 10)
	@test dL0.(xx) == -0.5 .* ones(length(xx))
	@test dL1.(xx) == 0.5 .* ones(length(xx))


	# constant functions
	L = constantFunctions()
	xx = LinRange(-1, 1, 10)
	@test L.(xx) == ones(length(xx))

	# basis function matrix
	N_matrix, dN_matrix = constructBasisFunctionMatrixLinearLagrange(evalPts = xx)
	R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts = xx)

	@test size(N_matrix) == size(dN_matrix) == (2, length(xx))
	@test R_matrix == ones(1, length(xx))
end


@testset "Gauss-Legendre quadrature rule" begin
	x, w = GaussLegendreQuadRule()

	@test x ≈ [-1 / sqrt(3); 1 / sqrt(3)]
	@test w ≈ [1; 1]

	x, w = GaussLegendreQuadRule(interval = [0, 1], numQuadPts = 3)

	@test sum(w) ≈ 1.0
	@test x ≈ [0.5 * (1 - sqrt(3 / 5)); 0.5; 0.5 * (1 + sqrt(3 / 5))]
end
