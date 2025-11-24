using Datasolver, Revise, LinearAlgebra, Test, Plots, ColorSchemes, LaTeXStrings, PGFPlotsX, JSON


paired_colors = colorschemes[:tableau_20]
single_colors = colorschemes[:tableau_10]
pgfplotsx()
default(size = (400, 300), markersize = 3, palette = single_colors, markerstrokewidth = 0.5)
# tickfont = font(16), guidefont = font(16), legendfont = font(18)
include("ANLP_solver.jl")

### simplified Kanno truss - 2-element truss

#region nonlinear strains and load-deflection with ANLP
A = 2000 / 1e6        # [m²]
bar_E = 1.654e+08   # [Pa]
βₛ = 1e-5           # scaling factor to improve the conditioning

α = 1.0

λ = 1
Fnodal = -400.0 * λ       # [N]
num_load_steps = 200

loadFac = LinRange(0.0, num_load_steps, num_load_steps + 1)

num_data_pts = 65

strain_limit = [3e-1;
	-3e-1]

dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1])


# geometry
node_vector = [
	[0, 0],
	[3.6, 0],
	[0, 3.6],
]

constrained_dofs = [
	(1, 1),
	(1, 2),
	(3, 1),
	(3, 2),
]

connections = [
	(1, 2),
	(3, 2),
]

# force
force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2

# truss problem
initProblem = TrussProblem(
	A,
	force,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
);

# ANLP results
resultsANLP = solveANLP(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	YoungModulus = bar_E,
	scaleFacYoungModulus = βₛ,
	NR_max_iter = 100,
	qrFactorized = false,
);


# plot load-deflection curve (deflection of node 2)
u2x = [resultsANLP.u[i][3] for i in 1:num_load_steps];
u2y = [resultsANLP.u[i][4] for i in 1:num_load_steps];


plot(-u2x, loadFac[2:end], label = "ux-node 2")
plot!(-u2y, loadFac[2:end], label = "uy-node 2")

plot!(framestyle = :box, xlabel = "uh", ylabel = "F[N]")

plot!(ylims = [0, 168])


# dataset
ii = 100;       # num_load_steps
scatter(dataset.E, dataset.S / βₛ, label = "dataset", framestyle = :box, xlabel = "strain", ylabel = "stress")

scatter!(resultsANLP.e[ii], resultsANLP.s[ii], marker = :rect, markersize = 4, label = "(eh,sh),ANLP")


# plot!(ylims = [-3e7, 3e7])


# Figure 6a:

# plot deformed structure at chosen load step
uh = resultsANLP.u[ii]
ux1 = uh[1:2:end]
uy1 = uh[2:2:end]

sc = 1.0
plot(0, 0, framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linecolor = :black)

	xN = [node_vector[i1][1] + ux1[i1] * sc, node_vector[i2][1] + ux1[i2] * sc]
	yN = [node_vector[i1][2] + uy1[i1] * sc, node_vector[i2][2] + uy1[i2] * sc]

	plot!(xN, yN, linecolor = :royalblue)
end

plot!(framestyle = :box, xlabel = "x", ylabel = "y", legend = false)

savefig("fig/kanno_trussSimp_deformed_structure.tex")

#endregion

# Figure 13 e/f

#region nonlinear strain + nonlinear dataset
A = 2000 / 1e6        # [m²]
βₛ = 1e-5           # scaling factor to improve the conditioning

α = 1.0

λ = 100
num_load_steps = 5


Fnodal = -400.0 * λ       # [N]

random_init_data = false
init_indices = nothing

num_data_pts = 65

strain_limit = [3e-1;
	-3e-1]
eSc = 10
smax = 3.5e7


stressFunc(x) = βₛ * smax .* (2 ./ (1 + exp(-eSc .* x)) - 1)
dataset = create_dataset(num_data_pts, stressFunc, strain_limit[2], strain_limit[1]);


# geometry
node_vector = [
	[0, 0],
	[3.6, 0],
	[0, 3.6],
]

constrained_dofs = [
	(1, 1),
	(1, 2),
	(3, 1),
	(3, 2),
]

connections = [
	(1, 2),
	(3, 2),
]

# force Vector
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2


# truss problem
initProblem = TrussProblem(
	A,
	force,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
);


# ADM
t1 = time();
resultsADM = Datasolver.directSolverNonLinearBarA(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	dataset = dataset,
	scaleFactorDataConst = βₛ,
	random_init_data = random_init_data,
	init_indices = init_indices,
	verbose = true,
	QRfactorized = false,
);
elapsed_timeADM = time() - t1

uh = resultsADM.u[end]
ux2 = uh[1:2:end]
uy2 = uh[2:2:end]

eh2 = resultsADM.e[end]
sh2 = resultsADM.s[end]


# GoADM
t2 = time();
resultsGoADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	dataset = dataset,
	scaleFactorDataConst = βₛ,
	random_init_data = random_init_data,
	init_indices = init_indices,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	verbose = true,
	QRfactorized = false,
);
elapsed_timeGoADM = time() - t2


uh = resultsGoADM.u[end]
ux3 = uh[1:2:end]
uy3 = uh[2:2:end]

eh3 = resultsGoADM.e[end]
sh3 = resultsGoADM.s[end]


#MINLP

if isfile("examples/kanno_trussSimp_minlp_solveresults.json")
	resultsMINLP = JSON.parsefile("examples/kanno_trussSimp_minlp_solveresults.json")

else
	initProblemMINLP = TrussProblem(
		A,
		force * βₛ,
		connections,
		α,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)
	t3 = time()
	resultsMINLP = Datasolver.NLP_solver(initProblemMINLP, dataset, use_L1_norm = false, use_data_bounds = true)
	elapsed_timeMINLP = time() - t3
	# save results to json
	open("examples/kanno_trussSimp_minlp_solveresults.json", "w") do io
		JSON.print(io, resultsMINLP)
	end
end




uh = resultsMINLP["u"][end]
ux4 = uh[1:2:end]
uy4 = uh[2:2:end]

eh4 = resultsMINLP["e"][end]
sh4 = resultsMINLP["s"][end] ./ βₛ


## plots
# dataset
# TODO find better markers
scatter(dataset.E, dataset.S / βₛ, label = "dataset", framestyle = :box, xlabel = "strain", ylabel = "stress")

# scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ, marker = :xcross, markersize = 5, markerstrokewidth = 2, label = "(etilde,stilde), ADM")

scatter!(resultsADM.e[end], resultsADM.s[end], marker = :circ, markersize = 4, label = L"$(e_h,s_h)$, ADM")

# scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ, marker = :cross, markersize = 5, markerstrokewidth = 2, label = "(etilde,stilde), GO-ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end], marker = :utriangle, markersize = 4, label = L"$(e_h,s_h)$, GO-ADM")

# scatter!(resultsMINLP.E[end], resultsMINLP.S[end] / βₛ, marker = :star5, markersize = 5, markerstrokewidth = 2, label = "(etilde,stilde), MINLP")

scatter!(resultsMINLP["e"][end], resultsMINLP["s"][end] / βₛ, marker = :diamond, markersize = 3, label = L"$(e_h,s_h)$, MINLP")


savefig("fig/kanno_trussSimp_dataset_nonlinE_nonlinData.tex")


# Figure 13a
# plot deformed structure
sc = 1.0

plot(0, 0, framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, label = i == 1 ? "Original" : nothing, linecolor = :black)

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, label = i == 1 ? "Deformed ADM" : nothing, linecolor = single_colors[3])

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, label = i == 1 ? "Deformed GO-ADM" : nothing, linecolor = single_colors[2])

	xN = [node_vector[i1][1] + ux4[i1] * sc, node_vector[i2][1] + ux4[i2] * sc]
	yN = [node_vector[i1][2] + uy4[i1] * sc, node_vector[i2][2] + uy4[i2] * sc]

	plot!(xN, yN, label = i == 1 ? "Deformed MINLP" : nothing, linecolor = single_colors[1], linestyle = :dash)
end
plot!(xlabel = "x", ylabel = "y", legend = true)

plot!(ylims = [-3, 4], yticks = [-3, -2, -1, 0, 1, 2, 3])
plot!(ylims = [-3, 4], yticks = [-2, 0, 2, 4])

savefig("fig/kanno_trussSimp_nonlinE_nonlinData_phih_100F.tex")



# plot stress
num_ele = length(connections)
plot(1:num_ele+1, [sh2[1]; sh2], linetype = :steppre, label = "ADM")
plot!(1:num_ele+1, [sh3[1]; sh3], linetype = :steppre, label = "GO-ADM")
plot!(1:num_ele+1, [sh4[1]; sh4], linetype = :steppre, label = "MINLP")


plot!(
	framestyle = :box,
	size = (800, 600),
	xticks = (1.5:1:11, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]),
	xlabel = "Element number",
	ylabel = "Axial stress",
	legend = :bottomright,
)

plot!(ylims = [-3e7, 3e7], yticks = [-2e7, 0, 2e7])

savefig("fig/kanno_trussSimp_nonlinE_nonlinData_sh_100F.tex")


#endregion




### TIMING


function print_iters()
	initProblem_linear = TrussProblem(
		A,
		force,
		connections,
		0.0,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)
	initProblem_nonlin = TrussProblem(
		A,
		force,
		connections,
		1.0,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)



	println("Linear ADM random init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	println("NonLinear ADM random init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)


	println("Linear ADM zero init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	println("NonLinear ADM zero init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	println("Linear ADM nullspace init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	println("NonLinear ADM nullspace init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)


	############## GOADM ################

	println("Linear GOADM random init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	println("NonLinear GOADM random init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)


	println("Linear GOADM zero init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	println("NonLinear GOADM zero init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	println("Linear GOADM nullspace init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	println("NonLinear GOADM nullspace init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)


end


function get_timings()
	initProblem_linear = TrussProblem(
		A,
		force,
		connections,
		0.0,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)
	initProblem_nonlin = TrussProblem(
		A,
		force,
		connections,
		1.0,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)
	lin_ADM_rand_times = zeros(100)
	nonlin_ADM_rand_times = zeros(100)
	lin_ADM_zero_times = zeros(100)
	nonlin_ADM_zero_times = zeros(100)
	lin_ADM_nullspace_times = zeros(100)
	nonlin_ADM_nullspace_times = zeros(100)

	lin_GOADM_rand_times = zeros(100)
	nonlin_GOADM_rand_times = zeros(100)
	lin_GOADM_zero_times = zeros(100)
	nonlin_GOADM_zero_times = zeros(100)
	lin_GOADM_nullspace_times = zeros(100)
	nonlin_GOADM_nullspace_times = zeros(100)


	println("Linear ADM random init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsADM = Datasolver.directSolverNonLinearBarA(
			initProblem = initProblem_linear,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = true,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		lin_ADM_rand_times[i] = elapsed_time
	end
	println("NonLinear ADM random init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsADM = Datasolver.directSolverNonLinearBarA(
			initProblem = initProblem_nonlin,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = true,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		nonlin_ADM_rand_times[i] = elapsed_time
	end

	println("Linear ADM zero init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsADM = Datasolver.directSolverNonLinearBarA(
			initProblem = initProblem_linear,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		lin_ADM_zero_times[i] = elapsed_time
	end
	println("NonLinear ADM zero init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsADM = Datasolver.directSolverNonLinearBarA(
			initProblem = initProblem_nonlin,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		nonlin_ADM_zero_times[i] = elapsed_time
	end
	println("Linear ADM nullspace init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsADM = Datasolver.directSolverNonLinearBarA(
			initProblem = initProblem_linear,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		lin_ADM_nullspace_times[i] = elapsed_time
	end
	println("NonLinear ADM nullspace init:")
	Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsADM = Datasolver.directSolverNonLinearBarA(
			initProblem = initProblem_nonlin,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		nonlin_ADM_nullspace_times[i] = elapsed_time
	end

	############## GOADM ################

	println("Linear GOADM random init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
			initProblem = initProblem_linear,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = true,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		lin_GOADM_rand_times[i] = elapsed_time
	end
	println("NonLinear GOADM random init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = true,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
			initProblem = initProblem_nonlin,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = true,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		nonlin_GOADM_rand_times[i] = elapsed_time
	end

	println("Linear GOADM zero init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
			initProblem = initProblem_linear,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		lin_GOADM_zero_times[i] = elapsed_time
	end
	println("NonLinear GOADM zero init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
			initProblem = initProblem_nonlin,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = ones(Int, length(connections)) * (num_data_pts ÷ 2 + 1),
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		nonlin_GOADM_zero_times[i] = elapsed_time
	end
	println("Linear GOADM nullspace init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_linear,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
			initProblem = initProblem_linear,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		lin_GOADM_nullspace_times[i] = elapsed_time
	end
	println("NonLinear GOADM nullspace init:")
	Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem_nonlin,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = false,
		init_indices = nothing,
		verbose = true,
		QRfactorized = false,
	)

	for i in 1:100
		t1 = time()
		resultsGOADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
			initProblem = initProblem_nonlin,
			constrained_dofs_global = constrained_dofs,
			externalForce = force,
			num_load_steps = num_load_steps,
			loadFac = Vector(loadFac),
			dataset = dataset,
			scaleFactorDataConst = βₛ,
			random_init_data = false,
			init_indices = nothing,
			verbose = false,
			QRfactorized = false,
		)
		elapsed_time = time() - t1
		nonlin_GOADM_nullspace_times[i] = elapsed_time
	end


	initProblemMINLP_linear = TrussProblem(
		A,
		force * βₛ,
		connections,
		0.0,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)
	initProblemMINLP_nonlinear = TrussProblem(
		A,
		force * βₛ,
		connections,
		1.0,
		constrained_dofs,
		node_vector = node_vector,
		num_quad_pts = 2,
	)


	lin_MINLP_times = zeros(100)

	for i in 1:100
		t3 = time()
		resultsMINLP = Datasolver.NLP_solver(initProblemMINLP_linear, dataset, use_L1_norm = false, use_data_bounds = true)
		elapsed_timeMINLP = time() - t3
		lin_MINLP_times[i] = elapsed_timeMINLP
	end

	nonlin_MINLP_times = zeros(100)

	for i in 1:100
		t3 = time()
		resultsMINLP = Datasolver.NLP_solver(initProblemMINLP_nonlinear, dataset, use_L1_norm = false, use_data_bounds = true)
		elapsed_timeMINLP = time() - t3
		nonlin_MINLP_times[i] = elapsed_timeMINLP
	end
	return Dict(
		"lin_ADM_zero_times" => lin_ADM_zero_times,
		"nonlin_ADM_zero_times" => nonlin_ADM_zero_times,
		"lin_ADM_rand_times" => lin_ADM_rand_times,
		"nonlin_ADM_rand_times" => nonlin_ADM_rand_times,
		"lin_ADM_nullspace_times" => lin_ADM_nullspace_times,
		"nonlin_ADM_nullspace_times" => nonlin_ADM_nullspace_times,
		"lin_GOADM_zero_times" => lin_GOADM_zero_times,
		"nonlin_GOADM_zero_times" => nonlin_GOADM_zero_times,
		"lin_GOADM_rand_times" => lin_GOADM_rand_times,
		"nonlin_GOADM_rand_times" => nonlin_GOADM_rand_times,
		"lin_GOADM_nullspace_times" => lin_GOADM_nullspace_times,
		"nonlin_GOADM_nullspace_times" => nonlin_GOADM_nullspace_times,
		"linMINLP_times" => lin_MINLP_times,
		"nonlinMINLP_times" => nonlin_MINLP_times)
end

if isfile("examples/simple_examples/simple_truss_times.json")
	println("Loading existing results...")
	all_times = open("examples/simple_examples/simple_truss_times.json", "r") do f
		JSON.parse(f)
	end
else
	all_times = get_timings()
	open("examples/simple_examples/simple_truss_times.json", "w") do f
		JSON.print(f, all_times)
	end
end

open("examples/simple_examples/simple_truss_times_summary.txt", "w") do f
	write(f, "Mean times:\n")
	lines = []
	for (key, value) in all_times
		push!(lines, "$key: $(round(mean(value),digits=5))\n")
	end
	write(f, join(sort(lines)))

	write(f, "Median times:\n")
	lines = []
	for (key, value) in all_times
		push!(lines, "$key: $(round(median(value),digits=5))\n")
	end
	write(f, join(sort(lines)))

end

if !isfile("examples/simple_examples/iter_prints.txt")
	touch("examples/simple_examples/iter_prints.txt")
	open("examples/simple_examples/iter_prints.txt", "w") do f
		redirect_stdout(f) do
			print_iters()
		end
	end

	# Remove lines that are simply "Skip this trial, already computed"
	open("examples/simple_examples/iter_prints.txt", "r") do f
		lines = readlines(f)
		lines = filter(x -> !contains(x, "Skip this trial, already computed"), lines)
		println(length(lines))
		open("examples/simple_examples/iter_prints.txt", "w") do f
			write(f, join(lines, "\n"))
		end
	end
end




# @show mean(ADM_times), mean(GOADM_times), mean(MINLP_times)
# scatter(zeros(100), ADM_times)
# scatter!(zeros(100) .+ 1, GOADM_times)
# scatter!(zeros(100) .+ 2, MINLP_times)
### End timing

# Figure 14

#region nonlinear strain + nonlinear dataset + init opt
A = 2000 / 1e6
βₛ = 1e-5

α = 1.0

λ = 100
num_load_steps = 5

Fnodal = -400.0 * λ

num_data_pts = 65
strain_limit = [3e-1;
	-3e-1]
eSc = 10
smax = 4e7

stressFunc(x) = βₛ * smax .* (2 ./ (1 + exp(-eSc .* x)) - 1)
dataset = create_dataset(num_data_pts, stressFunc, strain_limit[2], strain_limit[1]);


# geometry
node_vector = [
	[0, 0],
	[3.6, 0],
	[0, 3.6],
]

constrained_dofs = [
	(1, 1),
	(1, 2),
	(3, 1),
	(3, 2),
]

connections = [
	(1, 2),
	(3, 2),
]

# force Vector
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2


# truss problem
initProblem = TrussProblem(
	A,
	force,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
);

# running through 3 initialization options
tt = zeros(3, 2);
nriter = zeros(3, 2);
admiter = zeros(3, 2);
compcost = zeros(num_load_steps, 3, 2);

for i in 1:3
	if i == 1
		# stress-free
		init_indices = Int64.(33 .* ones(length(connections)))
		random_init_data = false
	elseif i == 2
		# random
		init_indices = nothing
		random_init_data = true
	else
		# nullspace
		init_indices = nothing
		random_init_data = false
	end

	# ADM
	t1 = time()
	resultsADM = Datasolver.directSolverNonLinearBarA(
		initProblem = initProblem,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = random_init_data,
		init_indices = init_indices,
		verbose = true,
		QRfactorized = false,
	)
	tt[i, 1] = time() - t1

	# GoADM
	t2 = time()
	resultsGoADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		dataset = dataset,
		scaleFactorDataConst = βₛ,
		random_init_data = random_init_data,
		init_indices = init_indices,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		verbose = true,
		QRfactorized = false,
	)
	tt[i, 2] = time() - t2

	# collect other metrics
	nriter[i, 1] = sum(sum.(resultsADM.NRiter))
	nriter[i, 2] = sum(sum.(resultsGoADM.NRiter))

	admiter[i, 1] = sum(resultsADM.ADMiter)
	admiter[i, 2] = sum(resultsGoADM.ADMiter)

	cc = 0
	for j in 1:num_load_steps
		cc += resultsADM.ADMiter[j]
		compcost[j, i, 1] = resultsADM.cost[cc]
	end
	compcost[:, i, 2] = resultsGoADM.cost
end






# MINLP
if !isfile("examples/kanno_trussSimp_load_step_MINLP_results.json")
	all_MINLP_results = []
	for load_step in 1:num_load_steps
		initProblemMINLP = TrussProblem(
			A,
			force * loadFac[load_step+1] * βₛ,
			connections,
			α,
			constrained_dofs,
			node_vector = node_vector,
			num_quad_pts = 2,
		)
		push!(all_MINLP_results, Datasolver.NLP_solver(initProblemMINLP, dataset, use_L1_norm = false, use_data_bounds = true))
	end
	open("examples/kanno_trussSimp_load_step_MINLP_results.json", "w") do f
		JSON.print(f, all_MINLP_results)
	end
else
	all_MINLP_results = open("examples/kanno_trussSimp_load_step_MINLP_results.json", "r") do f
		JSON.parse(f)
	end
end


# plot
lsty = [:solid, :dash, :dashdot]

plot(xlabel = "load step", ylabel = "value of the cost function", framestyle = :box)

for i in 1:3
	plot!(compcost[:, i, 1], label = "ADM, init opt $i", linestyle = lsty[i])
	plot!(compcost[:, i, 2], label = "GO-ADM, init opt $i", linestyle = lsty[i])
end


costs = [all_MINLP_results[i]["cost"] for i in 1:num_load_steps]


plot!(collect(1:num_load_steps),
	[all_MINLP_results[i]["cost"][1] for i in 1:num_load_steps],
	label = "MINLP", linestyle = :dot, linecolor = :black)

plot!(yscale = :log10)


plot!(ylims = (5e-6, 1e-3))

savefig("fig/kanno_trussSimp_costFunc_nonlinE_nonlinData.tex")


#endregion
