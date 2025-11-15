using Datasolver, Revise, LinearAlgebra, Test, Plots, JSON

include("ANLP_solver.jl")


#---------------------------------------------------------------------
#region linear strain
### benchmark the Kanno truss with linear solution computed with 
### https://valdivia.staff.jade-hs.de/fachwerk_en.html


# material
A = 2000 / 1e6        # [m²]
bar_E = 1.622e+09   # [Pa]
βₛ = 1e-4           # scaling factor of bar_E to improve the conditioning

α = 1.0         # 0: linear strain     1: nonlinear strain

# load steps and factors
Fnodal = -400.0       # [N]

num_load_steps = 1
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

# number data point and strain limit

strain_limit = [5e-4;
	-5e-4]

# geometry
node_vector = [
	[0.0, 0.0],
	[3.6, 0.0],
	[0.0, 3.6],
]

constrained_dofs = [
	(1, 1),
	(1, 2),
	(3, 1),
	(3, 2),
]

connections = [
	(1, 2),
	(2, 3),
]



#---------------------------------------------------------------------
#region nonlinear strains and comparing with ANLP

α = 1.0
βₛ = 1e-4
Fnodal = -400.0 * 500       # [N]


# force
force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2

# truss problem
initProblem = TrussProblem(
	A,
	force .* βₛ,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
)


# ADM results
begin
	num_data_pts = 65
	strain_limit = [2e-1;
		-2e-1]


	dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1], noise_magnitude = 0.1)


	# if isfile("simple_truss_results_GOADM.json")
	println("Running solver and storing results...")
	resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBar(initProblem, dataset)
	open("simple_truss_results_GOADM.json", "w") do f
		JSON.print(f, resultsGOADM)
	end
	println("Loading existing results...")
	resultsGOADM = open("simple_truss_results_GOADM.json", "r") do f
		JSON.parse(f)
	end
	# else

	# end


	# if isfile("simple_truss_results.json")
	println("Running solver and storing results...")
	resultsMINLP = Datasolver.NLP_solver(initProblem, dataset, use_L1_norm = false, use_data_bounds = true)
	@show resultsMINLP.u
	open("simple_truss_results.json", "w") do f
		JSON.print(f, resultsMINLP)
	end
	println("Loading existing results...")
	resultsMINLP = open("simple_truss_results.json", "r") do f
		JSON.parse(f)
	end
	# else

	# end

	# if isfile("simple_truss_results_ADM.json")
	println("Running solver and storing results...")
	resultsADM = Datasolver.directSolverNonLinearBar(initProblem, dataset)
	@show resultsADM.u
	open("simple_truss_results_ADM.json", "w") do f
		JSON.print(f, resultsADM)
	end

	println("Loading existing results...")
	resultsADM = open("simple_truss_results_ADM.json", "r") do f
		JSON.parse(f)
	end
	# else

	# end
	# resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBar(initProblem, dataset)

	resultsMINLP["s"] = resultsMINLP["s"] ./ βₛ
	resultsMINLP["S"] = resultsMINLP["S"] ./ βₛ
	resultsGOADM["S"] = resultsGOADM["S"] ./ βₛ
	resultsGOADM["s"] = resultsGOADM["s"] ./ βₛ
	resultsADM["S"] = resultsADM["S"] ./ βₛ
	resultsADM["s"] = resultsADM["s"] ./ βₛ
end

deformation_scale = 1

node_matrix = hcat(node_vector...)'
deformation_matrix = reshape(resultsMINLP["u"][1], 2, length(node_vector))'
deformation_matrix_goadm = reshape(resultsGOADM["u"][1], 2, length(node_vector))'
deformation_matrix_adm = reshape(resultsADM["u"][1], 2, length(node_vector))'
new_node_matrix = node_matrix + deformation_matrix * deformation_scale
new_node_matrix_goadm = node_matrix + deformation_matrix_goadm * deformation_scale
new_node_matrix_adm = node_matrix + deformation_matrix_adm * deformation_scale
scatter(node_matrix[:, 1], node_matrix[:, 2], label = "Original nodes", legend = :topright)
scatter!(new_node_matrix[:, 1], new_node_matrix[:, 2], label = "Deformed nodes", legend = :topright, marker = :xcross)
scatter!(new_node_matrix_goadm[:, 1], new_node_matrix_goadm[:, 2], label = "Deformed nodes GOADM", legend = :topright, marker = :diamond)
scatter!(new_node_matrix_adm[:, 1], new_node_matrix_adm[:, 2], label = "Deformed nodes ADM", legend = :topright, marker = :utriangle)

# scatter(dataset.E, dataset.S ./ βₛ, label = "Dataset", xlabel = "Strain", ylabel = "Stress [Pa]", legend = :topright)
# scatter!(resultsMINLP["e"], resultsMINLP["s"], label = "ADM solution", xlabel = "Strain", ylabel = "Stress [Pa]", legend = :topright)




scatter(dataset.E, dataset.S ./ βₛ, label = "Dataset", xlabel = "Strain", ylabel = "Stress [Pa]", legend = :topright)
scatter!(resultsGOADM["E"], resultsGOADM["S"], label = "GOADM solution", xlabel = "Strain", ylabel = "Stress [Pa]", legend = :topright)

scatter!(resultsMINLP["E"], resultsMINLP["S"], label = "MINLP solution", xlabel = "Strain", ylabel = "Stress [Pa]", legend = :topleft)

scatter!(resultsADM["E"], resultsADM["S"], label = "ADM solution", xlabel = "Strain", ylabel = "Stress [Pa]", legend = :topleft)
