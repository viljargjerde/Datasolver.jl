
# Test for get_final function
@testset "get_final Function Tests" begin
	solve_result = create_test_solve_results()
	push!(solve_result.E, [0.01, 0.02, 0.03])
	push!(solve_result.S, [100.0, 200.0, 300.0])

	final_values = get_final(solve_result)
	@test Dict(pairs(final_values)) == Dict(pairs((e = [], E = [0.01, 0.02, 0.03], s = [], S = [100.0, 200.0, 300.0], u = [], λ = [], μ = [], cost = [], equilibrium = [], compatibility = [], data_idx = Any[], solvetime = Any[])))
end

