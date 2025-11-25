
using CSV, DataFrames, XLSX, DelaunayTriangulation
using Datasolver, Revise, LinearAlgebra, Plots, LaTeXStrings, PGFPlotsX

pgfplotsx()

include("ANLP_solver.jl")


#region computing with linear and nonlinear strains
##### nonlinear (strain) 1D bar structure ######

bar_L = Float64(π)
A = π * 0.02^2          # [m²] 
bar_E = 7e10            # [Pa]
βₛ = 1e-6

α = 0.0         # 0: linear strain     1: nonlinear strain

ne = 8

numDataPts = 64+1

# manifactured solution and nonlinear force function
β = 0.15*π

uRef(x,β) = β * sin(π*x/bar_L)

eRef(x,α,β) = β * π / bar_L * cos(π*x/bar_L) * ( 1 + 0.5*α * β * π / bar_L * cos(π*x/bar_L) )

sRef(x,α,β) = bar_E * eRef(x,α,β)

force_func(α, β, x, λ) = λ .*
    [bar_E * A * ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_L))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_L)) / bar_L))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_L)) / bar_L) + 2) / bar_L^(2)); 0]

num_load_steps = 1

loadFac = LinRange(0,1.0,num_load_steps+1)

# mesh
h = bar_L/ne
node_vector = [ [(i-1)*h, 0] for i in 1:ne+1 ]

constrained_dofs = [
    (1, 1),
    (1, 2),
    (ne+1,1),
    (ne+1,2)
]

connections = [ (i, i+1) for i in 1:ne ]

# strain limit = max of eRef + safety increase
xx = 0:bar_L/1000:bar_L
eMax = maximum(abs(x[1]) for x in eRef.(xx,α,β))
strain_limit = 1.5 .* [eMax;
                       -eMax]

dataset = create_dataset(numDataPts, x -> bar_E*βₛ * x, strain_limit[2], strain_limit[1])


initProblem = TrussProblem(
    A,
    [0],
    connections,
    α,
    constrained_dofs,
    node_vector = node_vector,
    num_quad_pts = 2,
    force_func = force_func
)

result = Datasolver.directSolverNonLinearBarA(
        initProblem=initProblem,
        constrained_dofs_global=constrained_dofs,
        externalForce = (x,λ) -> force_func(α, β, x, λ),
        dataset=dataset,
        scaleFactorDataConst=βₛ,
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true
);
  
result2 = Datasolver.greedyLocalSearchSolverNonLinearBarA(
        initProblem=initProblem,
        constrained_dofs_global=constrained_dofs,
        externalForce = (x,λ) -> force_func(α, β, x, λ),
        dataset=dataset,
        scaleFactorDataConst=βₛ,
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true
);

checkThermomechanicalConsistency(results=result)

checkThermomechanicalConsistency(results=result2)

# taking results of the last ADM iter
uh = result.u[end]
uxh = uh[1:2:end]

eh = result.e[end]


uh2 = result2.u[end]
uxh2 = uh2[1:2:end]

eh2 = result2.e[end]


# plots

# force function 
forceX_func(α, β, x, λ) = λ *
     ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_L))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_L)) / bar_L))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_L)) / bar_L) + 2) / bar_L^(2));

plot(xx, forceX_func.(0.0, β, xx, 1.0), linewidth=2, linecolor=:black, label=L"$\alpha = 0.0$")
plot!(xx, forceX_func.(1.0, β, xx, 1.0), linewidth=2, linecolor=:royalblue, label=L"$\alpha = 1.0$")

plot!(legend=:topright, dpi=150, framestyle=:box, size=(800,600), xlabel=L"$\xi$ [m]", ylabel=L"$f (\xi)$ [N]", tickfont=font(24), guidefont=font(24), legendfont=font(24))

savefig("/scratch/ddcm/elsarticle/figs/1DBar_fx.pdf")


# dataset
scatter(dataset.E, dataset.S / βₛ / 1e9, marker=:circ, markercolor=:gray, markersize=5, markeralpha=0.7, markerstrokealpha=0, label=L"\tilde{y}", dpi=150, framestyle=:box, size=(800,600), xlabel="strain [-]", ylabel="stress [GPa]", tickfont=font(24), guidefont=font(24))

scatter!(result.E[end], result.S[end] / βₛ / 1e9, marker=:utriangle, markersize=10, markercolor=:crimson, markeralpha=0.3, markerstrokecolor=:crimson, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, ADM")
scatter!(result2.E[end], result2.S[end] / βₛ / 1e9, marker=:rect, markersize=10, markercolor=:forestgreen, markeralpha=0.3, markerstrokecolor=:forestgreen, markerstrokealpha=1, label=L"$\tilde{y}_h^*$, GO-ADM")

scatter!(result.e[end], result.s[end] / 1e9, marker=:cross, markersize=12, markercolor=:crimson, markeralpha=1, markerstrokecolor=:crimson, markerstrokealpha=1, label=L"$y_h$, ADM")
scatter!(result2.e[end], result2.s[end] / 1e9, marker=:cross, markersize=12, markercolor=:forestgreen, markeralpha=1, markerstrokecolor=:forestgreen, markerstrokealpha=1, label=L"$y_h$, GO-ADM")

plot!(legendfont=font(24), legend=:bottomright)



if α == 0
    plot!(xlims=(-0.75,0.75),xticks=[-0.8,-0.4,0,0.4,0.8])
    plot!(ylims=(-55,55))

    savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_dataset.pdf")
else
    plot!(xlims=(-0.9,0.9),xticks=[-0.8,-0.4,0,0.4,0.8])
    plot!(ylims=(-65,65))

    savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_dataset.pdf")
end



# ux
plot(xx, uRef.(xx,β), linewidth=2, linecolor=:black, label="Reference solution")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], uxh, linewidth=2, linecolor=:crimson,label="ADM")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], uxh2, linewidth=2, linecolor=:forestgreen,label="GO-ADM")

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel=L"$\xi$ [m]", ylabel=L"$u_{h,x}$", tickfont=font(24), guidefont=font(24), legendfont=font(24), legend=:bottom)


if α == 0
    savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_uh.pdf")
else
    savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_uh.pdf")
end


# axial strain
plot(xx,eRef.(xx,α,β), linewidth=2, linecolor=:black, label="Reference solution")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], [eh[1];eh], linewidth=2,linetype=:steppre, label="ADM", linecolor=:crimson)
plot!([node_vector[i][1] for i in 1:initProblem.num_node], [eh2[1];eh2], linewidth=2,linetype=:steppre,label="GO-ADM", linecolor=:forestgreen)

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel=L"$\xi$ [m]", ylabel=L"$e_h$ [-]", tickfont=font(24), guidefont=font(24), legendfont=font(24), legend=:topright)


if α == 0
    savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_eh.pdf")
else
    savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_eh.pdf")
end


#endregion



#----------------------------------------------------------------
#region CONVERGENCE STUDY

βₛ = 1e-8

# loop over num_ele and numDataPts
N_datapoints = [2^n for n in 4:10]
N_elements = [2^n for n in 2:8]

# allocation
l2e = zeros(length(N_datapoints),length(N_elements),2)
l2eGA = zeros(length(N_datapoints),length(N_elements),2)


for (aa, α) in enumerate([0.0,1.0])
    # ADM
    for (i,N_d) in enumerate(N_datapoints), (j,N_e) in enumerate(N_elements)
        @show N_d, N_e

        eMax = maximum(abs(x[1]) for x in eRef.(xx,α,β))
        strain_limit = 1.5 .* [eMax;
                            -eMax]

        local dataset = create_dataset(N_d, x -> bar_E*βₛ * x, strain_limit[2], strain_limit[1])

        # mesh
        h = bar_L/N_e
        node_vector = [ [(i-1)*h, 0] for i in 1:N_e+1 ]

        constrained_dofs = [
            (1, 1),
            (1, 2),
            (N_e+1,1),
            (N_e+1,2)
        ]

        connections = [ (i, i+1) for i in 1:N_e ]

        local initProblem = TrussProblem(A,
                                        [0],
                                        connections,
                                        α,
                                        constrained_dofs,
                                        node_vector = node_vector,
                                        num_quad_pts = 2,
                                        force_func = force_func)

        local result = Datasolver.directSolverNonLinearBarA(
                            initProblem=initProblem,
                            constrained_dofs_global=constrained_dofs,
                            externalForce = (x,λ) -> force_func(α, β, x, λ),
                            dataset=dataset,
                            scaleFactorDataConst=βₛ,
                            num_load_steps=num_load_steps,
                            loadFac=Vector(loadFac),
                            verbose=true,
                            NR_max_iter=20);

        l2e[i,j,aa] = Datasolver.relL2err1D(problem=initProblem, uNodal=result.u[end], uAfunction=x->uRef(x,β))
    end


    # GoADM
    for (i,N_d) in enumerate(N_datapoints), (j,N_e) in enumerate(N_elements)

        N_d = N_datapoints[i]
        N_e = N_elements[j]

        @show N_d, N_e

        eMax = maximum(abs(x[1]) for x in eRef.(xx,α,β))
        strain_limit = 1.5 .* [eMax;
                            -eMax]
        
        local dataset = create_dataset(N_d, x -> bar_E*βₛ * x, strain_limit[2], strain_limit[1])

        # mesh
        h = bar_L/N_e
        node_vector = [ [(i-1)*h, 0] for i in 1:N_e+1 ]

        constrained_dofs = [
            (1, 1),
            (1, 2),
            (N_e+1,1),
            (N_e+1,2)
        ]

        connections = [ (i, i+1) for i in 1:N_e ]

        local initProblem = TrussProblem(A,
                                        [0],
                                        connections,
                                        α,
                                        constrained_dofs,
                                        node_vector = node_vector,
                                        num_quad_pts = 2,
                                        force_func = force_func)

        local result2 = Datasolver.greedyLocalSearchSolverNonLinearBarA(
                        initProblem=initProblem,
                        constrained_dofs_global=constrained_dofs,
                        externalForce = (x,λ) -> force_func(α, β, x, λ),
                        dataset=dataset,
                        scaleFactorDataConst=βₛ,
                        num_load_steps=num_load_steps,
                        loadFac=Vector(loadFac),
                        verbose=true,
                        NR_max_iter=20);

        l2eGA[i,j,aa] = Datasolver.relL2err1D(problem=initProblem, uNodal=result2.u[end], uAfunction=x->uRef(x,β))
    end
end

### plots
clims = [(-6.3,-3), (-4.2,-3)]
for aa in 1:2
    # ADM
    #lvs = collect(-4:0.15:-2)
    #lvs_labels = map(x -> "$x", lvs)

    contour(
            N_elements, N_datapoints, log10.(l2e[:,:,aa]),
            ylabel = "Number of data points",
            xlabel = "Number of elements",
            #levels = lvs,
            #levels_label=lvs_labels,
            clabels=false, 
            color=:inferno,
            fill=false,
            colorbar=true,
            clim=clims[aa],
            colorbar_ticks = [-6, -5, -4, -3],
            colorbar_tickfontsize = 20,
            linewidth=3,
            #colorbar_title = L"Relative $L^2$ error $(\log_{10})$",
            #colorbar_titlefontsize = 24,
            #colorbar_orientation = :horizontal,
            scale = :log10,
            framestyle = :box,
            dpi=150, 
            size=(800,600), 
            tickfont=font(20), 
            guidefont=font(24),
            legendfont=font(20)
    )


    if aa == 1
        savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_convergenceADM.pdf")
    else
        savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_convergenceADM.pdf")
    end



    # GoADM
    #lvs = collect(-4:0.15:-2)
    #lvs_labels = map(x -> "$x", lvs)

    contour(
            N_elements, N_datapoints, log10.(l2eGA[:,:,aa]),
            ylabel = "Number of data points",
            xlabel = "Number of elements",
            #levels = lvs,
            #levels_label=lvs_labels,
            clabels=false, 
            color=:inferno,
            fill=false,
            colorbar=true,
            clim=clims[aa],
            colorbar_ticks = [-6, -5, -4, -3],
            colorbar_tickfontsize = 20,
            linewidth=3,
            #colorbar_title = L"Relative $L^2$ error $(\log_{10})$",
            #colorbar_titlefontsize = 24,
            #colorbar_orientation = :horizontal,
            scale = :log10,
            framestyle = :box,
            dpi=150, 
            size=(800,600), 
            tickfont=font(20), 
            guidefont=font(24),
            legendfont=font(20)
    )


    if aa == 1
        savefig("/scratch/ddcm/elsarticle/figs/1DlinBar_convergenceGoADM.pdf")
    else
        savefig("/scratch/ddcm/elsarticle/figs/1DnonlinBar_convergenceGoADM.pdf")
    end
end


#endregion
