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
