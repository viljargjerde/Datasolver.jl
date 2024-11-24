
function create_test_dataset()
	E = [0.01, 0.02, 0.03, 0.04]
	S = [100.0, 150.0, 200.0, 250.0]
	C = sum(S ./ (E .+ 1e-10)) / length(S)
	return Dataset(E, S, C)
end


# Test for the Dataset struct
@testset "Dataset Struct Tests" begin
	dataset = create_test_dataset()
	@test dataset.E == [0.01, 0.02, 0.03, 0.04]
	@test dataset.S == [100.0, 150.0, 200.0, 250.0]
	@test dataset.C ≈ sum(dataset.S ./ (dataset.E .+ 1e-10)) / length(dataset.S)
end



# Test for create_dataset function
@testset "Utils create_dataset min max stress Function Tests" begin
	N = 10
	strain_stress_relation = x -> 100.0 * x
	min_strain = -0.2
	max_strain = 0.1
	min_stress = 0.0
	max_stress = 10.0

	dataset = create_dataset(N, strain_stress_relation, min_strain, max_strain, min_stress, max_stress)

	@test length(dataset.E) == N
	@test length(dataset.S) == N
	@test dataset.C ≈ sum(dataset.S ./ (dataset.E .+ 1e-10)) / length(dataset.S)
	@test minimum(dataset.E) >= min_strain
	@test maximum(dataset.E) <= max_strain
	@test minimum(dataset.S) >= min_stress - 1e10
	@test maximum(dataset.S) <= max_stress + 1e-10
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
