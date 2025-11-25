
using CSV, DataFrames, XLSX, DelaunayTriangulation
using Datasolver, Revise, LinearAlgebra, Plots

include("ANLP_solver.jl")

# cross section
dd = 208/1000           # rope diameter [m]
A = (dd/2)^2 * π
bar_L = 17010 / 1000    # rope length [m]

βₛ = 1e-6
ne = 16
Fnodal = 1.0

α = 1.0

# mesh
h = bar_L/ne
node_vector = [ [(i-1)*h, 0] for i in 1:ne+1 ]

constrained_dofs = [
    (1, 1),
    (1, 2)
];

connections = [ (i, i+1) for i in 1:ne ];

# force Vector
force = zeros(2 * length(node_vector));
force[end-1] = Fnodal


# read dataset
myfile = "examples/simple_examples/cyclic_master.xlsx"
mysheet = "PR05"
myrange = "B4:G16516"

Adata = XLSX.readdata(myfile, mysheet, myrange);

setIds = [1:166, 167:263, 264:359]      # nD =  1:357 (complete 1st cycle)
                                        #       1:166 (1st loading path)
                                        #       167:263 (1st deloading path)
                                        #       264:359 (2nd loading path)
numDataset = length(setIds)

# plot the complete dataset and loading steps (for overview)
tt, ee, ff = Float64.(Adata[:,1]), Float64.(Adata[:,end]) ./ 100, Float64.(Adata[:,2]);
ss = 1000ff ./ A;

plot(tt[1:931] ./ 60, ff[1:931] / 1e3, linewidth=2, linecolor=:darkorange, xlabel=L"$t$ [min]", ylabel="Nodal axial force [MN]", dpi=150, framestyle=:box, size=(800,600), tickfont=font(24), guidefont=font(24), label=false, xlims=(-0.1,3.1))

savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_realData_loads.pdf")



plot(ee[setIds[3][end]+1:end], ss[setIds[3][end]+1:end] / 1e6, linewidth=1, linecolor=:gray, label="Provided dataset")
plot!(ee[setIds[1]], ss[setIds[1]] / 1e6, linewidth=2, linecolor=:royalblue, label="1st active dataset")
plot!(ee[setIds[2]], ss[setIds[2]] / 1e6, linewidth=2, linecolor=:crimson, label="2nd active dataset")
plot!(ee[setIds[3]], ss[setIds[3]] / 1e6, linewidth=2, linecolor=:forestgreen, label="3rd active dataset")

plot!(xlims=[0.028,0.046], ylims=[24,144], xlabel="strain [-]", ylabel="stress [MPa]", dpi=150, framestyle=:box, size=(800,600), tickfont=font(24), guidefont=font(24), legend=:topleft, legendfont=font(24))


savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_realData_datasetAll.pdf")



# define init problem
initProblem = TrussProblem(
    A,
    force,
    connections,
    α,
    constrained_dofs,
    node_vector = node_vector,
    num_quad_pts = 2
);

# allocation for initial solution guess and solution for the 2nd and 3rd dataset
ndofs = Datasolver.get_ndofs(initProblem)

uhA, ehA, shA = zeros(ndofs[1]), zeros(ndofs[2]), zeros(ndofs[3])
μhA, λhA = zeros(ndofs[4]), zeros(ndofs[5])

uhG, ehG, shG = zeros(ndofs[1]), zeros(ndofs[2]), zeros(ndofs[3])
μhG, λhG = zeros(ndofs[4]), zeros(ndofs[5])

data_idsA = nothing
data_idsG = nothing

# allocation for results at each load step
nnAloadsteps = setIds[end][end];
cc_load = 0

uxh = zeros(nnAloadsteps,2)
eh = zeros(nnAloadsteps,2)
sh = zeros(nnAloadsteps,2)
etilde = zeros(nnAloadsteps,2)
stilde = zeros(nnAloadsteps,2)
compcost = zeros(nnAloadsteps,2)

maxNRiter = zeros(numDataset,2)      # 3 datasets
maxADMiter = zeros(numDataset,2)     # 3 dataset

# dataset
Fdata, tdata = [Float64[] for i in 1:numDataset], [Float64[] for i in 1:numDataset]
Sdata, Edata = [Float64[] for i in 1:numDataset], [Float64[] for i in 1:numDataset]

for (cc_set, nD) in enumerate(setIds)   

    Fdata[cc_set] = Adata[nD,2] .* 1000             # [N]
    Sdata[cc_set] = Fdata[cc_set] ./ A              # [N/m^2]
    Edata[cc_set] = Adata[nD,end] ./ 100            # [-]
    tdata[cc_set] = Adata[nD,1]                     # [s]
   
    # creat (raw) dataset
    dataset = Dataset(Edata[cc_set], Sdata[cc_set] .* βₛ);

    # load steps/factors
    num_load_steps = length(Fdata[cc_set])
    loadFac = [0;Fdata[cc_set]]     # add zero force for counting reason

    # solving
    # get solution from the last computation of the previous dataset
    xInitADM = [uhA; ehA; shA .* βₛ; μhA .* βₛ; λhA]
    xInitGoADM = [uhG; ehG; shG .* βₛ; μhG .* βₛ; λhG]

    if cc_set > 1
        data_idsA = Datasolver.find_closest_idx(dataset.S, shA.* βₛ)
        data_idsG = Datasolver.find_closest_idx(dataset.S, shG.* βₛ)
    end

    
    resultsADM = Datasolver.directSolverNonLinearBarA(
                                initProblem=initProblem,
                                constrained_dofs_global=constrained_dofs,
                                externalForce=force,
                                num_load_steps=num_load_steps,
                                loadFac=loadFac,
                                dataset=dataset,
                                init_indices = data_idsA,
                                scaleFactorDataConst=βₛ,
                                solution_guess1stload=xInitADM,
                                DD_max_iter=50,
                                verbose=true,
                                QRfactorized=true
    );

    resultsGoADM = Datasolver.greedyLocalSearchSolverNonLinearBarA(
            initProblem=initProblem,
            constrained_dofs_global=constrained_dofs,
            externalForce = force,
            dataset=dataset,
            scaleFactorDataConst=βₛ,
            init_indices=data_idsG,
            num_load_steps=num_load_steps,
            loadFac=loadFac,
            solution_guess1stload=xInitGoADM,
            DD_max_iter=50,
            verbose=true,
            QRfactorized=true
    );

    # getting last converged results from the dataset
    uhA = resultsADM.u[end]
    ehA = resultsADM.e[end]
    shA = resultsADM.s[end]
    μhA = resultsADM.μ[end]
    λhA = resultsADM.λ[end]

    uhG = resultsGoADM.u[end]
    ehG = resultsGoADM.e[end]
    shG = resultsGoADM.s[end]
    μhG = resultsGoADM.μ[end]
    λhG = resultsGoADM.λ[end]

    # storing results per load step
    cc = 0
    for j in 1:num_load_steps
        cc += resultsADM.ADMiter[j]

        # axial tip displacement
        uxh[cc_load+j,1] = resultsADM.u[cc][end-1]
        uxh[cc_load+j,2] = resultsGoADM.u[j][end-1]

        # axial strain and stress
        eh[cc_load+j,1] = resultsADM.e[cc][end]
        eh[cc_load+j,2] = resultsGoADM.e[j][end]

        sh[cc_load+j,1] = resultsADM.s[cc][end]
        sh[cc_load+j,2] = resultsGoADM.s[j][end]

        # selected data point
        etilde[cc_load+j,1] = resultsADM.E[cc][end]
        etilde[cc_load+j,2] = resultsGoADM.E[j][end]
        
        stilde[cc_load+j,1] = resultsADM.S[cc][end]
        stilde[cc_load+j,2] = resultsGoADM.S[j][end]

        # computational cost
        compcost[cc_load+j,1] = resultsADM.cost[cc]
    end
    compcost[nD,2] = resultsGoADM.cost
    cc_load = nD[end]

    # store max NR iter and ADM iter per dataset
    maxNRiter[cc_set,1] = maximum(maximum.(resultsADM.NRiter))
    maxNRiter[cc_set,2] = maximum(maximum.(resultsGoADM.NRiter))

    maxADMiter[cc_set,1] = maximum(resultsADM.ADMiter)
    maxADMiter[cc_set,2] = maximum(resultsGoADM.ADMiter)
end


maxNRiter
maxADMiter


# collect all data
tD, eD, sD, fD = Float64[], Float64[], Float64[], Float64[]
for cc_set in 1:numDataset
    tD = [tD;collect(tdata[cc_set])]
    eD = [eD;collect(Edata[cc_set])]
    sD = [sD;collect(Sdata[cc_set])]
    fD = [fD;collect(Fdata[cc_set])]
end


# plot
# dataset
for cc_set in 1:numDataset
    nD = setIds[cc_set]
    plot(eD, sD/1e6, linewidth=2, linecolor=:gray, label=L"\tilde{y}")
    #scatter!(Edata[cc_set], Sdata[cc_set]/1e6,label="dataset no.$cc_set")

    scatter!(etilde[nD,1], stilde[nD,1]./ βₛ/1e6, marker=:utriangle, markersize=5, markercolor=:crimson, markeralpha=0.3, markerstrokecolor=:crimson, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, ADM")

    scatter!(etilde[nD,2], stilde[nD,2]./ βₛ/1e6, marker=:rect, markersize=5, markercolor=:forestgreen, markeralpha=0.3, markerstrokecolor=:forestgreen, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, GO-ADM")

    scatter!(eh[nD,1], sh[nD,1]/1e6, marker=:cross, markersize=8, markercolor=:crimson, markeralpha=1, markerstrokecolor=:crimson, markerstrokealpha=1, label=L"$y_h$, ADM")
    
    scatter!(eh[nD,2], sh[nD,2]/1e6, marker=:cross, markersize=8, markercolor=:forestgreen, markeralpha=1, markerstrokecolor=:forestgreen, markerstrokealpha=1, label=L"$y_h$, GO-ADM")

    plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="strain [-]", ylabel="stress [MPa]", tickfont=font(24), guidefont=font(24))

    plot!(xlims=[0.028,0.046], ylims=[24,144], legendfont=font(24), legend=:bottomright)

    if α == 1
        savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_realData_datasetNo$cc_set.pdf")
    else
        savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_realData_datasetNo$cc_set.pdf")
    end
end


# cost function
plot(compcost[:,1], linewidth = 2, linecolor = :crimson, label = "ADM")
plot!(compcost[:,2], linewidth = 2, linecolor = :forestgreen, label = "GO-ADM")

plot!(yscale = :log10, xlabel = "load step", ylabel = L"dist$_G(\cdot)$", dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(24), guidefont = font(24), legendfont = font(24), legend = :bottomright)
plot!(ylims = (5e-11, 1e-3))


if α == 1
    savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_realData_costFunc.pdf")
else
    savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_realData_costFunc.pdf")
end



# ux-F load-deflection curve
plot(uxh[:,1], fD/1e6, linewidth=2, linecolor=:crimson, label="ADM")
plot!(uxh[:,2], fD/1e6, linewidth=2, linecolor=:forestgreen,label="GO-ADM")

plot!(xlabel = L"$u_{h,x}$ [m]", ylabel = L"$f$ [MN]", legend=:topleft,  dpi = 150, framestyle = :box, size = (800, 600), tickfont = font(24), guidefont = font(24), legendfont = font(24))


if α == 1
    savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_realData_uFcurve.pdf")
else
    savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_realData_uFcurve.pdf")
end


if α == 1
    plot!(xlims=[0.7,0.74], ylims=[4e6,4.8e6], legend=false)
    savefig("fig/1DnonlinBar_realData_uFcurve_zoomin.png")
else
    plot!(xlims=[0.74,0.76], ylims=[4e6,4.8e6], legend=false)
    savefig("fig/1DlinBar_realData_uFcurve_zoomin.png")
end