
using CSV, DataFrames, XLSX, DelaunayTriangulation
using Datasolver, Revise, LinearAlgebra, Plots, LaTeXStrings, PGFPlotsX

pgfplotsx()

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

num_data_pts = 129

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
SEminlp = [-0.13125  -2.0000108811461397e7;
			0.225  2.8284353095228143e7 ];


scatter(dataset.E, dataset.S / βₛ / 1e6, marker=:circ, markercolor=:gray, markersize=5, markeralpha=0.7, markerstrokealpha=0, label=L"\tilde{y}", dpi=150, framestyle=:box, size=(800,600), xlabel="strain [-]", ylabel="stress [MPa]", tickfont=font(24), guidefont=font(24))

scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ / 1e6, marker=:utriangle, markersize=10, markercolor=:crimson, markeralpha=0.3, markerstrokecolor=:crimson, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, ADM")

scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ / 1e6, marker=:rect, markersize=10, markercolor=:forestgreen, markeralpha=0.3, markerstrokecolor=:forestgreen, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, GO-ADM")
scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ / 1e6, marker=:diamond, markersize=8, markercolor=:royalblue, markeralpha=0.3, markerstrokecolor=:royalblue, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, MINLP")

scatter!(resultsADM.e[end], resultsADM.s[end] / 1e6, marker=:cross, markersize=12, markercolor=:crimson, markeralpha=1, markerstrokecolor=:crimson, markerstrokealpha=1, label=L"$y_h$, ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end] / 1e6, marker=:cross, markersize=12, markercolor=:forestgreen, markeralpha=1, markerstrokecolor=:forestgreen, markerstrokealpha=1, label=L"$y_h$, GO-ADM")

scatter!(SEminlp[:,1], SEminlp[:,2] / 1e6, marker=:star, markersize=10, markercolor=:royalblue, markeralpha=1, markerstrokecolor=:royalblue, markerstrokealpha=1, label=L"$y_h$, MINLP")

plot!(legend = :bottomright, legendfont = font(24))

plot!(xlims=(-0.31,0.31), xticks=[-0.3,-0.15,0,0.15,0.3])
plot!(ylims=(-33,33), yticks=[-30,-15,0,15,30])


savefig("/scratch/ddcm/elsarticle/figs/kanno_trussSimp_dataset_nonlinE_nonlinData.pdf")



# plot deformed structure
sc = 1.0

plot(0, 0)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black, label = i == length(connections) ? "Reference configuration" : nothing)

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson, label = i == length(connections) ? "ADM" : nothing)

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :forestgreen, label = i == length(connections) ? "GO-ADM" : nothing)

	plot!(xN, yN, linewidth = 2, linecolor = :royalblue, label = i == length(connections) ? "MINLP" : nothing)

	# xN = [node_vector[i1][1] + ux4[i1] * sc, node_vector[i2][1] + ux4[i2] * sc]
	# yN = [node_vector[i1][2] + uy4[i1] * sc, node_vector[i2][2] + uy4[i2] * sc]

	# plot!(xN, yN, linewidth = 1, linecolor = :blue, linestyle = :dash, label = i == length(connections) ? "MINLP" : nothing)

	# xN = [node_vector[i1][1] + ux5[i1] * sc, node_vector[i2][1] + ux5[i2] * sc]
	# yN = [node_vector[i1][2] + uy5[i1] * sc, node_vector[i2][2] + uy5[i2] * sc]

	# plot!(xN, yN, linewidth = 1, linecolor = :green, linestyle = :dash, label = i == length(connections) ? "ADM2" : nothing)

	# xN = [node_vector[i1][1] + ux6[i1] * sc, node_vector[i2][1] + ux6[i2] * sc]
	# yN = [node_vector[i1][2] + uy6[i1] * sc, node_vector[i2][2] + uy6[i2] * sc]

	# plot!(xN, yN, linewidth = 1, linecolor = :orange, linestyle = :dot, label = i == length(connections) ? "GOADM2" : nothing)
end
plot!(xlims=(-0.1,4), ylims = (-2.5,4), legend = :topright, dpi = 150, framestyle = :box, size = (800, 600), xlabel = L"$x$ [m]", ylabel = L"$y$ [m]", tickfont = font(24), guidefont = font(24), legendfont = font(24))

plot!(yticks=[-2,0,2,4])

savefig("/scratch/ddcm/elsarticle/figs/kanno_trussSimp_nonlinE_nonlinData_phih_100F.pdf")



# plot stress
num_ele = length(connections)

# extract MINLP from tikz pic
sh4 = [-2.0000108811461397e7;
	   -2.0000108811461397e7;
		2.8284353095228143e7];

plot(1:num_ele+1, [sh2[1]; sh2] /1e6, linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:num_ele+1, [sh3[1]; sh3] /1e6, linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)
plot!(1:num_ele+1, sh4 /1e6, linewidth = 2, linetype = :steppre, label = "MINLP", linecolor = :royalblue)

plot!(
	dpi = 150,
	framestyle = :box,
	size = (800, 600),
	xticks = (1.5:1:3, ["1", "2"]),
	xlabel = "Element number",
	ylabel = "Axial stress [MPa]",
	tickfont = font(24),
	guidefont = font(24),
	legendfont = font(24),
	legend = :bottomright,
)

plot!(ylims = [-30, 30], yticks = [-20, 0, 20])

savefig("/scratch/ddcm/elsarticle/figs/kanno_trussSimp_nonlinE_nonlinData_sh_100F.pdf")


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



#region replot cost function using results from tikz

costFadm = [1.627384165217586e-5  1.627384165217586e-5	1.627384165217586e-5;
            0.00011290670765705271  0.00011290670765705271	0.00011290670765705271;
            0.0005610790376587109  0.0005610790376587109 0.0005610790376587109;
            0.00023914659178074735 0.00023914659178074735 0.00023914659178074735;
            0.0004135528832466316 0.0004135528832466316 0.0004135528832466316]


costFgoadm = [1.627384165217586e-5  1.627384165217586e-5 1.627384165217586e-5;
              0.00011290670765705271  0.00011290670765705271	0.00011290670765705271;
              0.00018099882809965854  0.00018099882809965854	0.00018099882809965854;
              5.179075034032132e-5    5.179075034032132e-5	5.179075034032132e-5;
              4.4058504682587046e-5  4.4058504682587046e-5	4.4058504682587046e-5]


costFminlp = [1.6268702412517996e-5  ;
              0.00011288604745991304 ;
              6.771459900875588e-5  ;
              7.493688061437686e-6  ;
              4.28930757950118e-5  ]



lsty = [:solid, :dash, :dashdotdot]
lbs = ["stress-free", "random", "structure-specific"]

plot(xlabel = "load step", ylabel = L"dist$_G(\cdot)$", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(24), guidefont = font(24))

for i in 1:3
	ll = string("ADM, ", lbs[i])
	plot!(costFadm[:, i], linewidth = 2, linecolor = :crimson, label = ll, linestyle = lsty[i])

	ll = string("GO-ADM, ", lbs[i])
	plot!(costFgoadm[:, i], linewidth = 2, linecolor = :forestgreen, label = ll, linestyle = lsty[i])

	if i==3
		plot!(costFminlp, linewidth = 2, linecolor = :royalblue, label = "MINLP")
	end
end
plot!(yscale = :log10)

plot!(ylims = (3e-9, 1e-3), legend = :bottom, legendfont = font(24))

plot!(yticks = [1e-8, 1e-6, 1e-4, 1e-2, 1e0])


savefig("/scratch/ddcm/elsarticle/figs/kanno_trussSimp_costFunc_nonlinE_nonlinData.pdf")


#endregion
