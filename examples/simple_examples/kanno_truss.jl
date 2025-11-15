
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

α = 0.0         # 0: linear strain     1: nonlinear strain

# load steps and factors
Fnodal = -400.0       # [N]

num_load_steps = 1
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

# number data point and strain limit
num_data_pts = 65

strain_limit = [5e-4;
	-5e-4]

# geometry
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
	force,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
)

dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1])


results = Datasolver.directSolverNonLinearBarA(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	dataset = dataset,
	scaleFactorDataConst = βₛ,
	verbose = true,
	QRfactorized = false,
);


# Go-ADM
results2 = Datasolver.greedyLocalSearchSolverNonLinearBarA(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	dataset = dataset,
	scaleFactorDataConst = βₛ,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	verbose = true,
	QRfactorized = false,
);

dataset_unscaled = create_dataset(num_data_pts, x -> bar_E * x, strain_limit[2], strain_limit[1])


results3 = Datasolver.NLP_solver(initProblem, dataset_unscaled, use_L1_norm = false, use_data_bounds = true, timelimit = 600.0)

checkThermomechanicalConsistency(results = results)
checkThermomechanicalConsistency(results = results2)
checkThermomechanicalConsistency(results = results3)


# extract solution
uh = results.u[end]
uxh = uh[1:2:end]
uyh = uh[2:2:end]

eh = results.e[end]
sh = results.s[end]

uh = results2.u[end]
uxh2 = uh[1:2:end]
uyh2 = uh[2:2:end]

eh2 = results2.e[end]
sh2 = results2.s[end]


uh = results3.u[end]
uxh3 = uh[1:2:end]
uyh3 = uh[2:2:end]

eh3 = results3.e[end]
sh3 = results3.s[end]

# linear solution from online truss caculator
uRef = [
	 0.000      0.000;
	-9.084e-4  -2.222e-3;
	-1.174e-3  -4.858e-3;
	 0.000      0.000;
	 8.672e-4  -2.065e-3;
	 1.045e-3  -4.680e-3
]

sRef = [-4.093e+5
	-1.198e+5
	-2.697e+5
	2.960e+5
	70979
	-1.135e+5
	1.694e+5
	80249
	3.907e+5
	80249
]

eRef = sRef ./ bar_E


# dataset
scatter(dataset.E, dataset.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(results.E[end], results.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde),ADM")
scatter!(results.e[end], results.s[end], marker = :circ, markersize = 8, label = "(eh,sh),ADM")

scatter!(results2.E[end], results2.S[end] / βₛ, marker = :cross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde),GO-ADM")
scatter!(results2.e[end], results2.s[end], marker = :utriangle, markersize = 8, label = "(eh,sh),GO-ADM")

scatter!(results3.E[end], results3.S[end], marker = :cross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde),MINLP")
scatter!(results3.e[end], results3.s[end], marker = :utriangle, markersize = 8, label = "(eh,sh),MINLP")


plot!(legendfont = font(18))

savefig("fig/kanno_truss_dataset_linE.png")


# plot deformed structure
sc = 5e1
plot(0, 0, dpi = 150, size = (800, 600), framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black)

	xN = [node_vector[i1][1] + uxh[i1] * sc, node_vector[i2][1] + uxh[i2] * sc]
	yN = [node_vector[i1][2] + uyh[i1] * sc, node_vector[i2][2] + uyh[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson)

	xN = [node_vector[i1][1] + uxh2[i1] * sc, node_vector[i2][1] + uxh2[i2] * sc]
	yN = [node_vector[i1][2] + uyh2[i1] * sc, node_vector[i2][2] + uyh2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :forestgreen)

	xN = [node_vector[i1][1] + uRef[i1, 1] * sc, node_vector[i2][1] + uRef[i2, 1] * sc]
	yN = [node_vector[i1][2] + uRef[i1, 2] * sc, node_vector[i2][2] + uRef[i2, 2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :royalblue)


	xN = [node_vector[i1][1] + uxh3[i1] * sc, node_vector[i2][1] + uxh3[i2] * sc]
	yN = [node_vector[i1][2] + uyh3[i1] * sc, node_vector[i2][2] + uyh3[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :green)

	xN = [node_vector[i1][1] + uRef[i1, 1] * sc, node_vector[i2][1] + uRef[i2, 1] * sc]
	yN = [node_vector[i1][2] + uRef[i1, 2] * sc, node_vector[i2][2] + uRef[i2, 2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :green)
end
plot!(legend = false, dpi = 150, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18))


savefig("fig/kanno_truss_linE_phih.png")


# plot stress
plot(1:11, [sh[1]; sh], linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:11, [sh2[1]; sh2], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)
plot!(1:11, [sh3[1]; sh3], linewidth = 2, linetype = :steppre, label = "MINLP", linecolor = :yellow)
plot!(1:11, [sRef[1]; sRef], linewidth = 2, linetype = :steppre, label = "linear reference solution", linecolor = :royalblue)

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


savefig("fig/kanno_truss_linE_sh.png")


#endregion
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#region nonlinear strains and comparing with ANLP

α = 1.0
βₛ = 1e-5
Fnodal = -400.0 * 1500       # [N]

num_load_steps = 200
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

# force
force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

# truss problem
initProblem = TrussProblem(
	A,
	force,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
)


initProblem_scaled_f = TrussProblem(
	A,
	force .* βₛ,
	connections,
	α,
	constrained_dofs,
	node_vector = node_vector,
	num_quad_pts = 2,
)

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

uh = resultsANLP.u[end]
ux1 = uh[1:2:end]
uy1 = uh[2:2:end]

eh1 = resultsANLP.e[end]
sh1 = resultsANLP.s[end]


# ADM results
num_data_pts = 11

strain_limit = [3e-1;
	-3e-1]


dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1])
dataset_unscaled = create_dataset(num_data_pts, x -> bar_E * x, strain_limit[2], strain_limit[1])


resultsADM = Datasolver.directSolverNonLinearBarA(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	dataset = dataset,
	scaleFactorDataConst = βₛ,
	verbose = true,
	QRfactorized = false,
);

uh = resultsADM.u[end]
ux2 = uh[1:2:end]
uy2 = uh[2:2:end]

eh2 = resultsADM.e[end]
sh2 = resultsADM.s[end]


# Go-ADM
resultsGoADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
	initProblem = initProblem,
	constrained_dofs_global = constrained_dofs,
	externalForce = force,
	dataset = dataset,
	scaleFactorDataConst = βₛ,
	num_load_steps = num_load_steps,
	loadFac = Vector(loadFac),
	verbose = true,
	QRfactorized = false,
);

uh = resultsGoADM.u[end]
ux3 = uh[1:2:end]
uy3 = uh[2:2:end]

eh3 = resultsGoADM.e[end]
sh3 = resultsGoADM.s[end]


# MINLP
resultsMINLP = Datasolver.NLP_solver(initProblem_scaled_f, dataset, use_L1_norm = true, use_data_bounds = true)
using Dates
println("Finished MINLP at ", Dates.now())
resultsMINLP2 = Datasolver.NLP_solver(initProblem_scaled_f, dataset, use_L1_norm = false, use_data_bounds = true)
println("Finished MINLP2 at ", Dates.now())




uh = resultsMINLP.u[end]
ux4 = uh[1:2:end]
uy4 = uh[2:2:end]

eh4 = resultsMINLP.e[end]
sh4 = resultsMINLP.s[end]


# dataset
scatter(dataset.E, dataset.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde),ADM")

scatter!(resultsADM.e[end], resultsADM.s[end], marker = :circ, markersize = 8, label = "(eh,sh),ADM")

scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde),GO-ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end], marker = :circ, markersize = 8, label = "(eh,sh),GO-ADM")

plot!(legendfont = font(18))


savefig("fig/kanno_truss_dataset_nonlinE.png")



# cost function

cc = 0
costADM = zeros(num_load_steps)
for i in 1:num_load_steps
	global cc += resultsADM.ADMiter[i]
	costADM[i] = resultsADM.cost[cc]
end

plot(costADM, linewidth = 2, linecolor = :crimson, label = "ADM")
plot!(resultsGoADM.cost, linewidth = 2, linecolor = :forestgreen, label = "GO-ADM")

plot!(yscale = :log10, xlabel = "load step", ylabel = "value of the cost function", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(16), guidefont = font(16), legendfont = font(18))

plot!(ylims = (1e-4, 1e-2), legend = :bottomright)

savefig("fig/kanno_truss_costFunc_nonlinE.png")


# plot deformed structure
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

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson)

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :forestgreen)
end
plot!(dpi = 150, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18), legend = false)


savefig("fig/kanno_truss_nonlinE_phih_1500F.png")



# plot stress
plot(1:11, [sh1[1]; sh1], linewidth = 2, linetype = :steppre, label = "ANLP", linecolor = :royalblue)
plot!(1:11, [sh2[1]; sh2], linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:11, [sh3[1]; sh3], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)

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


savefig("fig/kanno_truss_nonlinE_sh_1500F.png")

#endregion



#region nonlinear strain + nonlinear dataset
A = 2000 / 1e6        # [m²]
# bar_E = 1.622e+09   # [Pa]
βₛ = 1e-5           # scaling factor of bar_E to improve the conditioning

α = 1.0

if α == 0
	λ = 1
	num_load_steps = 5
else
	λ = 1500
	num_load_steps = 200
end

Fnodal = -400.0 * λ       # [N]    1500

random_init_data = false
init_indices = nothing  # Int64.(33 .* ones(10))         # nothing

num_data_pts = 65
if α == 0.0
	strain_limit = [3e-4;
		-3e-4]
	eSc = 3e4
	smax = 5e5
else
	strain_limit = [3e-1;
		-3e-1]
	eSc = 30
	smax = 5e8
end

stressFunc(x) = βₛ * smax .* (2 ./ (1 + exp(-eSc .* x)) - 1)
dataset = create_dataset(num_data_pts, stressFunc, strain_limit[2], strain_limit[1]);


# geometry
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

# force Vector
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

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


## plots
# dataset
scatter(dataset.E, dataset.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), ADM")

scatter!(resultsADM.e[end], resultsADM.s[end], marker = :circ, markersize = 8, label = "(eh,sh), ADM")

scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ, marker = :cross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), GO-ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end], marker = :utriangle, markersize = 8, label = "(eh,sh), GO-ADM")

plot!(legendfont = font(14))

savefig("fig/kanno_truss_dataset_linE_nonlinData.png")
# savefig("fig/kanno_truss_dataset_linE_nonlinData_randomInit.png")
# savefig("fig/kanno_truss_dataset_linE_nonlinData_zeroInit.png")

savefig("fig/kanno_truss_dataset_nonlinE_nonlinData.png")
# savefig("fig/kanno_truss_dataset_nonlinE_nonlinData_randomInit.png")
# savefig("fig/kanno_truss_dataset_nonlinE_nonlinData_zeroInit.png")



# cost function
cc = 0
costADM = zeros(num_load_steps)
for i in 1:num_load_steps
	cc += resultsADM.ADMiter[i]
	costADM[i] = resultsADM.cost[cc]
end

plot(costADM, linewidth = 2, linecolor = :crimson, label = "ADM")
plot!(resultsGoADM.cost, linewidth = 2, linecolor = :forestgreen, label = "GO-ADM")

plot!(yscale = :log10, xlabel = "load step", ylabel = "value of the cost function", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(16), guidefont = font(16), legendfont = font(18))

# plot!(ylims=(1e-8,1e-5), legend=:topleft)

# plot!(ylims=(1e-4,2e-1), legend=:bottom)


# savefig("fig/kanno_truss_costFunc_linE_nonlinData.png")
# # savefig("fig/kanno_truss_costFunc_linE_nonlinData_randomInit.png")
# # savefig("fig/kanno_truss_costFunc_linE_nonlinData_zeroInit.png")


# savefig("fig/kanno_truss_costFunc_nonlinE_nonlinData.png")
# # savefig("fig/kanno_truss_costFunc_nonlinE_nonlinData_randomInit.png")
# # savefig("fig/kanno_truss_costFunc_nonlinE_nonlinData_zeroInit.png")



# plot deformed structure
if α == 0
	sc = 50 * num_load_steps
else
	sc = 1.0
end
plot(0, 0, dpi = 150, size = (800, 600), framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black)

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson)

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :forestgreen)
end
plot!(dpi = 150, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18), legend = false)


savefig("fig/kanno_truss_linE_nonlinData_phih_F.png")
savefig("fig/kanno_truss_nonlinE_nonlinData_phih_1500F.png")



# plot stress
plot(1:11, [sh2[1]; sh2], linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:11, [sh3[1]; sh3], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)

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

savefig("fig/kanno_truss_linE_nonlinData_sh_F.png")
savefig("fig/kanno_truss_nonlinE_nonlinData_sh_1500F.png")



#endregion



#region nonlinear strain + nonlinear dataset + init opt
A = 2000 / 1e6
bar_E = 1.622e+09
βₛ = 1e-5

α = 0.0

if α == 0
	λ = 1
	num_load_steps = 5
else
	λ = 1500
	num_load_steps = 200
end

Fnodal = -400.0 * λ

num_data_pts = 65
if α == 0.0
	strain_limit = [3e-4;
		-3e-4]
	eSc = 3e4
	smax = 5e5
else
	strain_limit = [3e-1;
		-3e-1]
	eSc = 30
	smax = 5e8
end

stressFunc(x) = βₛ * smax .* (2 ./ (1 + exp(-eSc .* x)) - 1)
dataset = create_dataset(num_data_pts, stressFunc, strain_limit[2], strain_limit[1]);


# geometry
node_vector = [
	[0, 0],
	[3.6, 0],
	[2 * 3.6, 0],
	[0, 3.6],
	[3.6, 3.6],
	[2 * 3.6, 3.6],
];

constrained_dofs = [
	(1, 1),
	(1, 2),
	(4, 1),
	(4, 2),
];

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
];

# force Vector
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

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
		init_indices = Int64.(33 .* ones(10))
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
	plot!(compcost[:, i, 1], linewidth = 2, linecolor = :crimson, label = "ADM, init opt $i", linestyle = lsty[i])
	plot!(resultsGoADM.cost, linewidth = 2, linecolor = :forestgreen, label = "GO-ADM, init opt $i", linestyle = lsty[i])
end

plot!(yscale = :log10)

plot!(ylims = (1e-8, 1e-5), legend = :topleft)

plot!(ylims = (3e-4, 2e-1), legend = :bottom)


savefig("fig/kanno_truss_costFunc_linE_nonlinData.png")

savefig("fig/kanno_truss_costFunc_nonlinE_nonlinData.png")


#endregion



#region nonlinear strain + unsymmetric dataset
A = 2000 / 1e6        # [m²]
bar_E = 1.622e+09   # [Pa]
βₛ = 1e-5           # scaling factor of bar_E to improve the conditioning

α = 1.0

if α == 0
	λ = 1
	num_load_steps = 5
else
	λ = 1500
	num_load_steps = 200
end

Fnodal = -400.0 * λ

num_data_pts = 65
if α == 0.0
	strain_limit = [3e-4;
		-3e-4]
	eSc = 3e4
	smax = 5e5
else
	strain_limit = [3e-1;
		-3e-1]
	eSc = 30
	smax = 5e8
end

stressFunc(x) = βₛ * smax .* (2 ./ (1 + exp(-eSc .* x)) - 1)

Neplus = Int64(ceil(0.8 * num_data_pts))
Neminus = Int64(num_data_pts - Neplus) + 1

dataE = [collect(range(strain_limit[2], 0.0, Neminus))[1:end-1];
	collect(range(0.0, strain_limit[1], Neplus))]

dataset = Dataset(dataE, stressFunc.(dataE))

# geometry
node_vector = [
	[0, 0],
	[3.6, 0],
	[2 * 3.6, 0],
	[0, 3.6],
	[3.6, 3.6],
	[2 * 3.6, 3.6],
];

constrained_dofs = [
	(1, 1),
	(1, 2),
	(4, 1),
	(4, 2),
];

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
];

# force Vector
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

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

# initialization options
tt = zeros(3, 2);
nriter = zeros(3, 2);
admiter = zeros(3, 2);
compcost = zeros(num_load_steps, 3, 2);

i = 2       # 1: stress-free    2: random   3: nullspace

if i == 1
	# stress-free
	init_indices = Int64.(33 .* ones(10))
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
tt[i, 1] = time() - t1

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
tt[i, 2] = time() - t2

uh = resultsGoADM.u[end]
ux3 = uh[1:2:end]
uy3 = uh[2:2:end]

eh3 = resultsGoADM.e[end]
sh3 = resultsGoADM.s[end]


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


## plots
# dataset
scatter(dataset.E, dataset.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), ADM")

scatter!(resultsADM.e[end], resultsADM.s[end], marker = :circ, markersize = 8, label = "(eh,sh), ADM")

scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ, marker = :cross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), GO-ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end], marker = :utriangle, markersize = 8, label = "(eh,sh), GO-ADM")

plot!(legendfont = font(14))

savefig("fig/kanno_truss_dataset_nonlinE_unsymData.png")


# cost function
lsty = [:solid, :dash, :dashdot]

plot(xlabel = "load step", ylabel = "value of the cost function", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(16), guidefont = font(16), legendfont = font(18))

for i in 1:3
	plot!(compcost[:, i, 1], linewidth = 2, linecolor = :crimson, label = "ADM, init opt $i", linestyle = lsty[i])
	plot!(resultsGoADM.cost, linewidth = 2, linecolor = :forestgreen, label = "GO-ADM, init opt $i", linestyle = lsty[i])
end
plot!(yscale = :log10)

plot!(ylims = (1e-4, 7e-1), legend = :bottomright)
plot!(yticks = [1e-4, 1e-3, 1e-2, 1e-1])

savefig("fig/kanno_truss_costFunc_nonlinE_unsymData.png")



# plot deformed structure
if α == 0
	sc = 50 * num_load_steps
else
	sc = 1.0
end
plot(0, 0, dpi = 150, size = (800, 600), framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black)

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson)

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :forestgreen)
end
plot!(dpi = 150, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18), legend = false)

savefig("fig/kanno_truss_nonlinE_unsymData_phih_1500F.png")



# plot stress
plot(1:11, [sh2[1]; sh2], linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:11, [sh3[1]; sh3], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)

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

savefig("fig/kanno_truss_nonlinE_unsymData_sh_1500F.png")


#endregion



#region nonlinear strain + noisy dataset
A = 2000 / 1e6        # [m²]
bar_E = 1.622e+09   # [Pa]
βₛ = 1e-5           # scaling factor of bar_E to improve the conditioning

α = 1.0

if α == 0
	λ = 1
	num_load_steps = 5
else
	λ = 1500
	num_load_steps = 200
end

Fnodal = -400.0 * λ

num_data_pts = 87
if α == 0.0
	strain_limit = [3e-4;
		-3e-4]
	eSc = 3e4
	smax = 5e5
else
	strain_limit = [3e-1;
		-3e-1]
	eSc = 30
	smax = 5e8
end

stressFunc(x) = βₛ * smax .* (2 ./ (1 + exp(-eSc .* x)) - 1)

dataset = create_dataset(num_data_pts, stressFunc, strain_limit[2], strain_limit[1]);

datasetNoisy = addGaussNoise(dataset = dataset, standardDev = 0.07, add2Sonly = false, add2Eonly = false, add2both = true)

ids = checkDatasetThermomechanicalConsistency(dataset = datasetNoisy)

# remove noise in the inconsistent data points
datasetNoisy.E[ids] = dataset.E[ids]
datasetNoisy.S[ids] = dataset.S[ids]


# geometry
node_vector = [
	[0, 0],
	[3.6, 0],
	[2 * 3.6, 0],
	[0, 3.6],
	[3.6, 3.6],
	[2 * 3.6, 3.6],
];

constrained_dofs = [
	(1, 1),
	(1, 2),
	(4, 1),
	(4, 2),
];

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
];

# force Vector
loadFac = LinRange(0.0, 1.0, num_load_steps + 1)

force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

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

# initialization options
tt = zeros(3, 2);
nriter = zeros(3, 2);
admiter = zeros(3, 2);
compcost = zeros(num_load_steps, 3, 2);

for i ∈ 1:3       # 1: stress-free    2: random   3: nullspace

	if i == 1
		# stress-free
		init_indices = Int64.(33 .* ones(10))
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
		dataset = datasetNoisy,
		scaleFactorDataConst = βₛ,
		random_init_data = random_init_data,
		init_indices = init_indices,
		verbose = true,
		QRfactorized = false,
	)
	tt[i, 1] = time() - t1

	uh = resultsADM.u[end]
	ux2 = uh[1:2:end]
	uy2 = uh[2:2:end]

	eh2 = resultsADM.e[end]
	sh2 = resultsADM.s[end]


	# GoADM
	t2 = time()
	resultsGoADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
		initProblem = initProblem,
		constrained_dofs_global = constrained_dofs,
		externalForce = force,
		dataset = datasetNoisy,
		scaleFactorDataConst = βₛ,
		random_init_data = random_init_data,
		init_indices = init_indices,
		num_load_steps = num_load_steps,
		loadFac = Vector(loadFac),
		search_iters = 200,
		verbose = true,
		QRfactorized = false,
	)
	tt[i, 2] = time() - t2

	uh = resultsGoADM.u[end]
	ux3 = uh[1:2:end]
	uy3 = uh[2:2:end]

	eh3 = resultsGoADM.e[end]
	sh3 = resultsGoADM.s[end]


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


## plots
# dataset
scatter(datasetNoisy.E, datasetNoisy.S / βₛ, label = "dataset", dpi = 150, framestyle = :box, size = (800, 600), xlabel = "strain", ylabel = "stress", tickfont = font(16), guidefont = font(16))

scatter!(resultsADM.E[end], resultsADM.S[end] / βₛ, marker = :xcross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), ADM")

scatter!(resultsADM.e[end], resultsADM.s[end], marker = :circ, markersize = 8, label = "(eh,sh), ADM")

scatter!(resultsGoADM.E[end], resultsGoADM.S[end] / βₛ, marker = :cross, markersize = 10, markerstrokewidth = 2, label = "(etilde,stilde), GO-ADM")

scatter!(resultsGoADM.e[end], resultsGoADM.s[end], marker = :utriangle, markersize = 8, label = "(eh,sh), GO-ADM")

plot!(legendfont = font(12))

savefig("fig/kanno_truss_dataset_nonlinE_noisyData.png")



# plot deformed structure
if α == 0
	sc = 50 * num_load_steps
else
	sc = 1.0
end
plot(0, 0, dpi = 150, size = (800, 600), framestyle = :box)

for i in 1:length(connections)
	i1, i2 = connections[i]
	xN = [node_vector[i1][1], node_vector[i2][1]]
	yN = [node_vector[i1][2], node_vector[i2][2]]

	plot!(xN, yN, linewidth = 2, linecolor = :black)

	xN = [node_vector[i1][1] + ux2[i1] * sc, node_vector[i2][1] + ux2[i2] * sc]
	yN = [node_vector[i1][2] + uy2[i1] * sc, node_vector[i2][2] + uy2[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :crimson)

	xN = [node_vector[i1][1] + ux3[i1] * sc, node_vector[i2][1] + ux3[i2] * sc]
	yN = [node_vector[i1][2] + uy3[i1] * sc, node_vector[i2][2] + uy3[i2] * sc]

	plot!(xN, yN, linewidth = 2, linecolor = :forestgreen)
end
plot!(dpi = 150, framestyle = :box, size = (800, 600), xlabel = "x", ylabel = "y", tickfont = font(16), guidefont = font(16), legendfont = font(18), legend = false)


savefig("fig/kanno_truss_nonlinE_noisyData_phih_1500F.png")



# plot stress
plot(1:11, [sh2[1]; sh2], linewidth = 2, linetype = :steppre, label = "ADM", linecolor = :crimson)
plot!(1:11, [sh3[1]; sh3], linewidth = 2, linetype = :steppre, label = "GO-ADM", linecolor = :forestgreen)

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

savefig("fig/kanno_truss_nonlinE_noisyData_sh_1500F.png")



# cost function
lsty = [:solid, :dash, :dashdot]

plot(xlabel = "load step", ylabel = "value of the cost function", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(16), guidefont = font(16), legendfont = font(18))

for i in 1:3
	plot!(compcost[:, i, 1], linewidth = 2, linecolor = :crimson, label = "ADM, init opt $i", linestyle = lsty[i])
	plot!(resultsGoADM.cost, linewidth = 2, linecolor = :forestgreen, label = "GO-ADM, init opt $i", linestyle = lsty[i])
end
plot!(yscale = :log10)

plot!(ylims = (1e-5, 4e1), legend = :bottom)
plot!(yticks = [1e-5, 1e-3, 1e-1, 1e1])

savefig("fig/kanno_truss_costFunc_nonlinE_noisyData.png")


#endregion
