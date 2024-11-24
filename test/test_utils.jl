function create_test_solve_results()
	Φ = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
	return SolveResults(N_datapoints = 4, Φ = Φ)
end


# Test for the SolveResults struct
@testset "SolveResults Struct Tests" begin
	solve_result = create_test_solve_results()
	@test solve_result.N_datapoints == 4
	@test solve_result.Φ == [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
	@test length(solve_result.E) == 0 # Testing default initialization
end

# Test for get_rel_diff function
@testset "Utils get_rel_diff Function Tests" begin
	xs = [0.0, 0.5, 1.0]
	u_solved = [0.0, 0.25, 0.5]
	u_analytical = [0.0, 0.24, 0.48]

	rel_diff_1 = get_rel_diff(xs, u_solved, u_analytical)

	@test rel_diff_1 >= 0.0

	u_solved = [0.0, 0.25, 0.5]
	u_analytical = u_solved

	rel_diff = get_rel_diff(xs, u_solved, u_analytical)
	@test rel_diff == 0.0

	u_analytical = [0.0, 0.22, 0.46]
	rel_diff = get_rel_diff(xs, u_solved, u_analytical)
	@test rel_diff_1 < rel_diff

end

@testset "Utils extract_matrix_block Tests" begin
	row_sizes = [2, 4, 4]
	col_sizes = [3, 1, 5]
	A_11 = rand(2, 3)
	A_12 = rand(2, 1)
	A_13 = rand(2, 5)

	A_21 = rand(4, 3)
	A_22 = rand(4, 1)
	A_23 = rand(4, 5)

	A_31 = rand(4, 3)
	A_32 = rand(4, 1)
	A_33 = rand(4, 5)

	A = [
		A_11 A_12 A_13
		A_21 A_22 A_23
		A_31 A_32 A_33
	]

	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 1, 1) == A_11
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 2, 1) == A_21
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 3, 1) == A_31
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 1, 2) == A_12
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 2, 2) == A_22
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 3, 2) == A_32
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 1, 3) == A_13
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 2, 3) == A_23
	@test Datasolver.extract_matrix_block(A, row_sizes, col_sizes, 3, 3) == A_33

end
