using Datasolver, Revise, LinearAlgebra, Test, Plots

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
	[3.6, 3.6],
	[0.0, 3.6],
]

constrained_dofs = [
	(1, 1),
	(1, 2),
	(4, 1),
	(4, 2),
]

connections = [
	(1, 2),
	(2, 3),
	(2, 4),
	(3, 4),
	(1, 3),
]


# force
force = zeros(2 * length(node_vector))
force[4] = 2 * Fnodal   # [N]   - downward force at node 2

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

#---------------------------------------------------------------------

#---------------------------------------------------------------------
#region nonlinear strains and comparing with ANLP

α = 1.0
βₛ = 1e-5
Fnodal = -400.0 * 1500       # [N]


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
	num_data_pts = 129
	strain_limit = [3e-1;
		-3e-1]


	dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1])

	resultsMINLP = Datasolver.NLP_solver(initProblem, dataset, use_L1_norm = false, use_data_bounds = true)


end


using Dates
println("Finished MINLP at ", Dates.now())
# resultsMINLP2 = Datasolver.NLP_solver(initProblem, dataset_unscaled, use_L1_norm = false, use_data_bounds = true)
println("Finished MINLP2 at ", Dates.now())






### full truss:


node_vector = [
	[0, 0],
	[3.6, 0],
	[2 * 3.6, 0],
	[0, 3.6],
	[3.6, 3.6],
	[2 * 3.6, 3.6],
]

constrained_dofs = [
	(1, 1),
	(1, 2),
	(4, 1),
	(4, 2),
]

connections = [
	(1, 2),
	(2, 3),
	(1, 5),
	(2, 4),
	(2, 5),
	(2, 6),
	(3, 5),
	(3, 6),
	(4, 5),
	(5, 6),
]


# force
force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

# truss problem
initProblem = TrussProblem(
	A,
	force .* βₛ ./ 5,
	connections,
	0.0,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
)

begin
	num_data_pts = 129
	strain_limit = [3e-1;
		-3e-1]


	dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1])

	resultsMINLP = Datasolver.NLP_solver(initProblem, dataset, use_L1_norm = false, use_data_bounds = true)
	resultsADM = Datasolver.directSolverNonLinearBar(initProblem, dataset)
	# resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBar(initProblem, dataset)

end



resultsMINLP.s ./ βₛ

scatter(resultsMINLP.E[end], resultsMINLP.S[end] ./ βₛ, marker = :circle, markersize = 6, label = "(etilde,stilde),MINLP", linecolor = :orange)
scatter!(resultsMINLP.e[end], resultsMINLP.s[end] ./ βₛ, marker = :diamond, markersize = 6, label = "(eh,sh),MINLP", linecolor = :red)
