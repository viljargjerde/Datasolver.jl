
using Datasolver, Revise, LinearAlgebra, Test, Plots

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


plot(-u2x, loadFac[2:end], linewidth = 2, linecolor = :royalblue, label = "ux-node 2")
plot!(-u2y, loadFac[2:end], linewidth = 2, linecolor = :crimson, label = "uy-node 2")

plot!(dpi = 150, framestyle = :box, size = (800, 600), xlabel = "uh", ylabel = "F[N]", tickfont = font(16), guidefont = font(16), legendfont = font(18))

plot!(ylims = [0, 168])


# dataset
ii = 100;       # num_load_steps
scatter(dataset.E, dataset.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(resultsANLP.e[ii], resultsANLP.s[ii], marker = :rect, markersize = 8, label = "(eh,sh),ANLP")

plot!(legendfont = font(18))
plot!(ylims = [-3e7, 3e7])



# plot deformed structure at chosen load step
uh = resultsANLP.u[ii]
ux1 = uh[1:2:end]
uy1 = uh[2:2:end]

sc = 1.0
plot(0, 0, dpi = 150, size = (800, 600), framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black)

	xN = [node_vector[i1][1] + ux1[i1] * sc, node_vector[i2][1] + ux1[i2] * sc]
	yN = [node_vector[i1][2] + uy1[i1] * sc, node_vector[i2][2] + uy1[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :royalblue)
end

plot!(dpi = 150, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18), legend = false)



#endregion



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


initProblem_scaled_f = TrussProblem(
	A,
	force .* βₛ,
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

t2 = time();
resultsADM2 = Datasolver.directSolverNonLinearBar(
	initProblem_scaled_f,
	dataset,
);
elapsed_timeADM = time() - t2

t2 = time();
resultsGOADM2 = Datasolver.greedyLocalSearchSolverNonLinearBar(
	initProblem_scaled_f,
	dataset,
);
elapsed_timeGOADM2 = time() - t2

uh = resultsADM.u[end]
ux2 = uh[1:2:end]
uy2 = uh[2:2:end]

eh2 = resultsADM.e[end]
sh2 = resultsADM.s[end]


uh = resultsADM2.u[end]
ux5 = uh[1:2:end]
uy5 = uh[2:2:end]

eh5 = resultsADM2.e[end]
sh5 = resultsADM2.s[end]

uh = resultsGOADM2.u[end]
ux6 = uh[1:2:end]
uy6 = uh[2:2:end]

eh6 = resultsGOADM2.e[end]
sh6 = resultsGOADM2.s[end]


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

resultsMINLP = Datasolver.NLP_solver(initProblem_scaled_f, dataset, use_L1_norm = false)

uh = resultsGoADM.u[end]
ux3 = uh[1:2:end]
uy3 = uh[2:2:end]

eh3 = resultsGoADM.e[end]
sh3 = resultsGoADM.s[end]


uh = resultsMINLP.u[end]
ux4 = uh[1:2:end]
uy4 = uh[2:2:end]

eh4 = resultsMINLP.e[end]
sh4 = resultsMINLP.s[end] ./ βₛ



## plots
# dataset
scatter(dataset.E, dataset.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), ADM")

scatter!(resultsADM.e[end], resultsADM.s[end], marker = :circ, markersize = 8, label = "(eh,sh), ADM")

scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ, marker = :cross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), GO-ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end], marker = :utriangle, markersize = 8, label = "(eh,sh), GO-ADM")


scatter!(resultsMINLP.E[end], resultsMINLP.S[end] ./ βₛ, marker = :diamond, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), MINLP")

scatter!(resultsMINLP.e[end], resultsMINLP.s[end] ./ βₛ, marker = :utriangle, markersize = 8, label = "(eh,sh), MINLP")

plot!(legendfont = font(14))

savefig("fig/kanno_trussSimp_dataset_nonlinE_nonlinData.png")



# plot deformed structure
sc = 1.0

plot(0, 0, dpi = 150, size = (800, 600), framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black, label = i == length(connections) ? "Original" : nothing)

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson, label = i == length(connections) ? "ADM" : nothing)

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, linewidth = 3, linecolor = :forestgreen, label = i == length(connections) ? "GO-ADM" : nothing)

	xN = [node_vector[i1][1] + ux4[i1] * sc, node_vector[i2][1] + ux4[i2] * sc]
	yN = [node_vector[i1][2] + uy4[i1] * sc, node_vector[i2][2] + uy4[i2] * sc]

	plot!(xN, yN, linewidth = 1, linecolor = :blue, linestyle = :dash, label = i == length(connections) ? "MINLP" : nothing)

	# xN = [node_vector[i1][1] + ux5[i1] * sc, node_vector[i2][1] + ux5[i2] * sc]
	# yN = [node_vector[i1][2] + uy5[i1] * sc, node_vector[i2][2] + uy5[i2] * sc]

	# plot!(xN, yN, linewidth = 1, linecolor = :green, linestyle = :dash, label = i == length(connections) ? "ADM2" : nothing)

	# xN = [node_vector[i1][1] + ux6[i1] * sc, node_vector[i2][1] + ux6[i2] * sc]
	# yN = [node_vector[i1][2] + uy6[i1] * sc, node_vector[i2][2] + uy6[i2] * sc]

	# plot!(xN, yN, linewidth = 1, linecolor = :orange, linestyle = :dot, label = i == length(connections) ? "GOADM2" : nothing)
end
plot!(dpi = 300, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18), legend = false)

plot!(ylims = [-3, 4], yticks = [-3, -2, -1, 0, 1, 2, 3])
plot!(ylims = [-3, 4], yticks = [-2, 0, 2, 4])
# plot!(ylims = [-3, 4], yticks = [-2, 0, 2, 4], legend = :outertopleft)

savefig("fig/kanno_trussSimp_nonlinE_nonlinData_phih_100F.png")



# plot stress
num_ele = length(connections)
plot(1:num_ele+1, [sh2[1]; sh2], linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:num_ele+1, [sh3[1]; sh3], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)
plot!(1:num_ele+1, [sh4[1]; sh4], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :blue)


plot!(
	dpi = 150,
	framestyle = :box,
	size = (800, 600),
	xticks = (1.5:1:11, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]),
	xlabel = "Element number",
	ylabel = "Axial stress",
	tickfont = font(16),
	guidefont = font(16),
	legendfont = font(18),
	legend = :bottomright,
)

plot!(ylims = [-3e7, 3e7], yticks = [-2e7, 0, 2e7])

savefig("fig/kanno_trussSimp_nonlinE_nonlinData_sh_100F.png")


#endregion



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


# plot
lsty = [:solid, :dash, :dashdot]

plot(xlabel = "load step", ylabel = "value of the cost function", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(16), guidefont = font(16), legendfont = font(18))

for i in 1:3
	plot!(compcost[:, i, 1], linewidth = 2, label = "ADM, init opt $i", linestyle = lsty[i])
	plot!(compcost[:, i, 2], linewidth = 2, label = "GO-ADM, init opt $i", linestyle = lsty[i])
end

plot!(yscale = :log10)


plot!(ylims = (5e-6, 1e-3), legend = :bottom)

savefig("fig/kanno_trussSimp_costFunc_nonlinE_nonlinData.png")


#endregion
