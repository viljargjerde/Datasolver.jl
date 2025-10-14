
using Datasolver, Revise, LinearAlgebra, Test, Plots


#---------------------------------------------------------------------
### benchmark the Kanno truss with linear solution computed with 
### https://valdivia.staff.jade-hs.de/fachwerk_en.html


# material
A = 2000/1e6        # [m²]
bar_E = 1.622e+09   # [Pa]
βₛ = 1e-4           # scaling factor of bar_E to improve the conditioning

α = 1.0         # 0: linear strain     1: nonlinear strain

# load steps and factors
Fnodal = -400.0       # [N]

num_load_steps = 1
loadFac = LinRange(0.0,1.0,num_load_steps+1)

# number data point and strain limit
num_data_pts = 65

strain_limit = [ 5e-4;
                -5e-4]

# geometry
node_vector = [
    [0,     0],
    [3.6,   0],
    [2*3.6, 0],
    [0,     3.6],
    [3.6,   3.6],
    [2*3.6, 3.6]
]  

constrained_dofs = [
    (1, 1),
    (1, 2),
    (4, 1),
    (4, 2)
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
    (5, 6)
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
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    dataset=dataset,
    scaleFactorDataConst=βₛ,
    verbose=true,
    QRfactorized=false
);


# Go-ADM
results = Datasolver.greedyLocalSearchSolverNonLinearBarA(
        initProblem=initProblem,
        constrained_dofs_global=constrained_dofs,
        externalForce = force,
        dataset=dataset,
        scaleFactorDataConst=βₛ,
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true,
        QRfactorized=false
);


checkThermomechanicalConsistency(results=results)

# extract solution
uh = results.u[end]
uxh = uh[1:2:end]
uyh = uh[2:2:end]

eh = results.e[end]
sh = results.s[end]


# linear solution from online truss caculator
uRef = [
 	 0.000 	    0.000;
 	-9.084e-4 	-2.222e-3;
 	-1.174e-3 	-4.858e-3;
 	 0.000 	    0.000;
 	 8.672e-4 	-2.065e-3;
 	 1.045e-3 	-4.680e-3
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

# plot deformed structure
sc = 5e1
plot(0,0,dpi=150,size=(800,600),framestyle=:box)

for i in 1:length(connections)
    i1,i2 = connections[i]
    xN = [node_vector[i1][1],node_vector[i2][1]]
    yN = [node_vector[i1][2],node_vector[i2][2]]

    plot!(xN,yN,linewidth=2,linecolor=:black)

    xN = [node_vector[i1][1]+uxh[i1]*sc,node_vector[i2][1]+uxh[i2]*sc]
    yN = [node_vector[i1][2]+uyh[i1]*sc,node_vector[i2][2]+uyh[i2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:crimson)

    xN = [node_vector[i1][1]+uRef[i1,1]*sc,node_vector[i2][1]+uRef[i2,1]*sc]
    yN = [node_vector[i1][2]+uRef[i1,2]*sc,node_vector[i2][2]+uRef[i2,2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:royalblue)
end
plot!(legend=false,dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="y", tickfont=font(16), guidefont=font(16),legendfont=font(18))


savefig("fig/kanno_truss_phih10F.png")


# plot stress
plot(1:11, [sh[1];sh], linewidth=2,linetype=:steppre,label="ddcm",linecolor=:crimson)
plot!(1:11, [sRef[1];sRef], linewidth=2,linetype=:steppre, label="lin. ref. sol.",linecolor=:royalblue)

plot!(dpi=150, framestyle=:box, size=(800,600), xticks=(1.5:1:11,["1","2","3","4","5","6","7","8","9","10"]), xlabel="Element number", ylabel="Axial stress", tickfont=font(16), guidefont=font(16),legendfont=font(18), legend=:bottomright)


savefig("fig/kanno_truss_sh10F.png")

#---------------------------------------------------------------------

#---------------------------------------------------------------------
### comparing with ANLP

include("ANLP_solver.jl")

βₛ = 1e-5 
Fnodal = -400.0*500       # [N]

num_load_steps = 1000
loadFac = LinRange(0.0,1.0,num_load_steps+1)

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

# ANLP results
resultsANLP = solveANLP(
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    YoungModulus = bar_E,
    scaleFacYoungModulus=βₛ,
    NR_max_iter = 100,
    qrFactorized = false
);

uh = resultsANLP.u[end]
ux1 = uh[1:2:end]
uy1 = uh[2:2:end]

eh1 = resultsANLP.e[end]
sh1 = resultsANLP.s[end]


# ADM results
num_data_pts = 129

strain_limit = [ 2e-1;
                -2e-1]


dataset = create_dataset(num_data_pts, x -> bar_E * βₛ * x, strain_limit[2], strain_limit[1])

resultsADM = Datasolver.directSolverNonLinearBarA(
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    dataset=dataset,
    scaleFactorDataConst=βₛ,
    verbose=true,
    QRfactorized=false
);

uh = resultsADM.u[end]
ux2 = uh[1:2:end]
uy2 = uh[2:2:end]

eh2 = resultsADM.e[end]
sh2 = resultsADM.s[end]


# Go-ADM
resultsGoADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
        initProblem=initProblem,
        constrained_dofs_global=constrained_dofs,
        externalForce = force,
        dataset=dataset,
        scaleFactorDataConst=βₛ,
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true,
        QRfactorized=false
);

uh = resultsGoADM.u[end]
ux3 = uh[1:2:end]
uy3 = uh[2:2:end]

eh3 = resultsGoADM.e[end]
sh3 = resultsGoADM.s[end]



# plot deformed structure
sc = 1.0
plot(0,0,dpi=150,size=(800,600),framestyle=:box)

for i in 1:length(connections)
    i1,i2 = connections[i]
    xN = [node_vector[i1][1],node_vector[i2][1]]
    yN = [node_vector[i1][2],node_vector[i2][2]]

    plot!(xN,yN,linewidth=2,linecolor=:black)

    xN = [node_vector[i1][1]+ux1[i1]*sc,node_vector[i2][1]+ux1[i2]*sc]
    yN = [node_vector[i1][2]+uy1[i1]*sc,node_vector[i2][2]+uy1[i2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:crimson)

    xN = [node_vector[i1][1]+ux2[i1]*sc,node_vector[i2][1]+ux2[i2]*sc]
    yN = [node_vector[i1][2]+uy2[i1]*sc,node_vector[i2][2]+uy2[i2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:royalblue)

    xN = [node_vector[i1][1]+ux3[i1]*sc,node_vector[i2][1]+ux3[i2]*sc]
    yN = [node_vector[i1][2]+uy3[i1]*sc,node_vector[i2][2]+uy3[i2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:forestgreen)
end
plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="y", tickfont=font(16), guidefont=font(16),legendfont=font(18), legend=false)


# plot stress
plot(1:11, [sh1[1];sh1], linewidth=2,linetype=:steppre,label="ANLP",linecolor=:crimson)
plot!(1:11, [sh2[1];sh2], linewidth=2,linetype=:steppre, label="ADM",linecolor=:royalblue)
plot!(1:11, [sh3[1];sh3], linewidth=2,linetype=:steppre, label="GoADM",linecolor=:forestgreen)

plot!(dpi=150, framestyle=:box, size=(800,600), xticks=(1.5:1:11,["1","2","3","4","5","6","7","8","9","10"]), xlabel="Element number", ylabel="Axial stress", tickfont=font(16), guidefont=font(16),legendfont=font(18), legend=:bottomright)
