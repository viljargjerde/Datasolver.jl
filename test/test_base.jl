using Test, SparseArrays, LinearAlgebra




@testset "LagrangePolynomials" begin
	# linear Lagrange polynomials
	L0, L1 = Datasolver.linearLagrangePolynomials([0.0, 1.0])

	@test L0(0) == 1.0
	@test L0(1) == 0

	@test L1(0) == 0.0
	@test L1(1) == 1


	# 1st derivative of linear Lagrange polynomials
	dL0, dL1 = Datasolver.compute1stDeriv4linearLagrangePolynomials([0.0, 1.0])
	xx = [0.25, 0.75]
	@test dL0.(xx) == -1 .* ones(length(xx))
	@test dL1.(xx) == 1 .* ones(length(xx))




	# basis function matrix
	N_mats, dN_mats = Datasolver.constructBasisFunctionMatrixLinearLagrange(1, xx)
	N_matrix = N_mats[1]
	dN_matrix = dN_mats[1]
	@test size(N_matrix) == size(dN_matrix) == (1, length(xx))
end


@testset "Gauss-Legendre quadrature rule" begin
	x, w = Datasolver.GaussLegendreQuadRule()

	@test x ≈ [-1 / sqrt(3); 1 / sqrt(3)]
	@test w ≈ [1; 1]

	x, w = Datasolver.GaussLegendreQuadRule(interval = [0, 1], numQuadPts = 3)

	@test sum(w) ≈ 1.0
	@test x ≈ [0.5 * (1 - sqrt(3 / 5)); 0.5; 0.5 * (1 + sqrt(3 / 5))]
end
