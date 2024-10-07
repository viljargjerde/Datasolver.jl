using LinearAlgebra
using Datasolver
using Test


# Defining some helper functions
function create_test_dataset()
	E = [0.01, 0.02, 0.03, 0.04]
	S = [100.0, 150.0, 200.0, 250.0]
	C = sum(S ./ (E .+ 1e-10)) / length(S)
	return Dataset(E, S, C)
end

function create_test_solve_results()
	Φ = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
	return SolveResults(N_datapoints = 4, Φ = Φ)
end


@testset "Datasolver.jl" begin
	# Write your tests here.
	# Test for the Dataset struct
	@testset "Dataset Struct Tests" begin
		dataset = create_test_dataset()
		@test dataset.E == [0.01, 0.02, 0.03, 0.04]
		@test dataset.S == [100.0, 150.0, 200.0, 250.0]
		@test dataset.C ≈ sum(dataset.S ./ (dataset.E .+ 1e-10)) / length(dataset.S)
	end


	# Test for the SolveResults struct
	@testset "SolveResults Struct Tests" begin
		solve_result = create_test_solve_results()
		@test solve_result.N_datapoints == 4
		@test solve_result.Φ == [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
		@test length(solve_result.E) == 0 # Testing default initialization
	end

	# Test for get_final function
	@testset "get_final Function Tests" begin
		solve_result = create_test_solve_results()
		push!(solve_result.E, [0.01, 0.02, 0.03])
		push!(solve_result.S, [100.0, 200.0, 300.0])

		final_values = get_final(solve_result)
		@test final_values == (e = [], E = [0.01, 0.02, 0.03], s = [], S = [100.0, 200.0, 300.0], u = [], η = [], μ = [], cost = [], balance = [], compatability = [])
	end

	# Test for remove_dofs(B::Matrix, node_coordinate_pairs) -> Matrix
	@testset "remove_dofs Matrix Tests" begin
		B = [1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0; 9.0 10.0 11.0 12.0]
		node_coordinate_pairs = [(2, 1)]

		new_B = Datasolver.remove_dofs(B, node_coordinate_pairs)
		@test new_B == [1.0 2.0 4.0; 5.0 6.0 8.0; 9.0 10.0 12.0]
	end

	# Test for remove_dofs(f::Vector, node_coordinate_pairs) -> Vector
	@testset "remove_dofs Vector Tests" begin
		f = [10.0, 20.0, 30.0, 40.0]
		node_coordinate_pairs = [(2, 1)]
		new_f = Datasolver.remove_dofs(f, node_coordinate_pairs)
		@test new_f == [10.0, 20.0, 40.0]

		node_coordinate_pairs = [(2, 2)]
		new_f = Datasolver.remove_dofs(f, node_coordinate_pairs)
		@test new_f == [10.0, 20.0, 30.0]

		node_coordinate_pairs = [(1, 1), (1, 2), (2, 2)]
		new_f = Datasolver.remove_dofs(f, node_coordinate_pairs)
		@test new_f == [30.0]
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


	# Test for create_dataset function
	@testset "Utils create_dataset min max stress Function Tests" begin
		N = 10
		strain_stress_relation = x -> 100.0 * x
		min_strain = -0.2
		max_strain = 0.1
		min_stress = 0.0
		max_stress = 10.0
		noise_magnitude = 0.05

		dataset = create_dataset(N, strain_stress_relation, min_strain, max_strain, min_stress, max_stress; noise_magnitude = noise_magnitude)

		@test length(dataset.E) == N
		@test length(dataset.S) == N
		@test dataset.C ≈ sum(dataset.S ./ (dataset.E .+ 1e-10)) / length(dataset.S)
		@test minimum(dataset.E) >= min_strain
		@test maximum(dataset.E) <= max_strain
		@test minimum(dataset.S) >= min_stress
		@test maximum(dataset.S) <= max_stress
	end

	@testset "Utils create_dataset Function Tests" begin
		N = 10
		strain_stress_relation = x -> 100.0 * x
		min_strain = -0.2
		max_strain = 0.1
		noise_magnitude = 0.05

		dataset = create_dataset(N, strain_stress_relation, min_strain, max_strain; noise_magnitude = noise_magnitude)

		@test length(dataset.E) == N
		@test length(dataset.S) == N
		@test dataset.C ≈ sum(dataset.S ./ (dataset.E .+ 1e-10)) / length(dataset.S)
		@test minimum(dataset.E) >= min_strain
		@test maximum(dataset.E) <= max_strain
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

end


