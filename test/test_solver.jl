

# Test for get_final function
@testset "get_final Function Tests" begin
	solve_result = create_test_solve_results()
	push!(solve_result.E, [0.01, 0.02, 0.03])
	push!(solve_result.S, [100.0, 200.0, 300.0])

	final_values = get_final(solve_result)
	@test final_values == (e = [], E = [0.01, 0.02, 0.03], s = [], S = [100.0, 200.0, 300.0], u = [], λ = [], μ = [], cost = [], balance = [], compatibility = [])
end



# Test for choose_closest_to function
@testset "choose_closest_to Function Tests" begin
	dataset = create_test_dataset()
	target_e = [0.02, 0.03]
	target_s = [150.0, 200.0]
	w = [1.0, 1.0]
	(e_values, s_values), cost = Datasolver.choose_closest_to(target_e, target_s, w, dataset)

	# Check choose_closest_to correctly finds the closest match
	@test e_values == [0.02, 0.03]
	@test s_values == [150.0, 200.0]
	@test cost >= 0.0 # Cost should be non-negative
end



# Test for datasolve function
@testset "datasolve Function Tests" begin
	connections = [[1, 2], [2, 3]]
	Φ = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
	A = [1.0, 1.0]
	dataset = create_test_dataset()
	f_with_dofs = [10.0, 20.0, 10.0, 20.0, 10.0, 20.0]
	f_without = [10.0, 20.0]
	fixed_dofs = [(1, 1), (1, 2), (3, 1), (3, 2)]

	result = datasolve(connections, Φ, A, dataset, f_with_dofs, fixed_dofs; initialization = true, max_iterations = 1000, verbose = false)
	result_without_dof = datasolve(connections, Φ, A, dataset, f_without, fixed_dofs; initialization = true, max_iterations = 1000, verbose = false)
	@test result.u[end] ≈ result_without_dof.u[end]
	@test typeof(result) == SolveResults # Test that it found a solution

end
