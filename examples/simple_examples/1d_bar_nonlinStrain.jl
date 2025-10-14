
using Datasolver, Revise, LinearAlgebra, Test, Plots, LaTeXStrings


##### nonlinear (strain) 1D bar structure ######

bar_L = Float64(π)
A = π * 0.02^2          # [m²] 
bar_E = 7e10            # [Pa]
sc = 1e-6

α = 1.0         # 0: linear strain     1: nonlinear strain

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

dataset = create_dataset(numDataPts, x -> bar_E*sc * x, strain_limit[2], strain_limit[1])


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

nonlin_result = Datasolver.directSolverNonLinearBarA(
        initProblem=initProblem,
        constrained_dofs_global=constrained_dofs,
        externalForce = (x,λ) -> force_func(α, β, x, λ),
        dataset=dataset,
        scaleFactorDataConst=sc,
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true
);
  
nonlin_result2 = Datasolver.greedyLocalSearchSolverNonLinearBarA(
        initProblem=initProblem,
        constrained_dofs_global=constrained_dofs,
        externalForce = (x,λ) -> force_func(α, β, x, λ),
        dataset=dataset,
        scaleFactorDataConst=sc,
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true
);

checkThermomechanicalConsistency(results=nonlin_result)

checkThermomechanicalConsistency(results=nonlin_result2)

# taking results of the last ADM iter
uh = nonlin_result.u[end]
uxh = uh[1:2:end]

eh = nonlin_result.e[end]


uh2 = nonlin_result2.u[end]
uxh2 = uh2[1:2:end]

eh2 = nonlin_result2.e[end]



# plots

# force function 
forceX_func(α, β, x, λ) = λ *
     ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_L))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_L)) / bar_L))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_L)) / bar_L) + 2) / bar_L^(2));

plot(xx, forceX_func.(0.0, β, xx, 1.0), linewidth=2, linecolor=:black,label="alpha = 0.0")
plot!(xx, forceX_func.(1.0, β, xx, 1.0), linewidth=2, linecolor=:magenta,label="alpha = 1.0")

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="fx/c [N]", tickfont=font(16), guidefont=font(16),legendfont=font(18))

savefig("fig/1DBar_fx.png")


# dataset
scatter(dataset.E, dataset.S / 1e4, label=false, dpi=150, framestyle=:box, size=(800,600), xlabel="strain", ylabel="stress x 1e-10", tickfont=font(16), guidefont=font(16))

savefig("fig/1DBar_dataset.png")


# ux
plot(xx, uRef.(xx,β), linewidth=2, linecolor=:black,label="uRef")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], uxh, linewidth=2, linecolor=:royalblue,label="uh-ADM")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], uxh2, linewidth=2, linecolor=:crimson,label="uh-GoADM")

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="Displacement ux", tickfont=font(16), guidefont=font(16),legendfont=font(18))
plot!(legend=:bottom)


savefig("fig/1DnonlinBar_uh.png")
savefig("fig/1DlinBar_uh.png")


# axial strain
plot(xx,eRef.(xx,α,β), linewidth=2, linecolor=:black,label="eRef")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], [eh[1];eh], linewidth=2,linetype=:steppre,label="eh-ADM", linecolor=:royalblue)
plot!([node_vector[i][1] for i in 1:initProblem.num_node], [eh2[1];eh2], linewidth=2,linetype=:steppre,label="eh-GoADM", linecolor=:crimson)

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="Axial strain", tickfont=font(16), guidefont=font(16),legendfont=font(18))


savefig("fig/1DnonlinBar_eh.png")
savefig("fig/1DlinBar_eh.png")



#----------------------------------------------------------------
# CONVERGENCE STUDY

# loop over num_ele and numDataPts

N_datapoints = [2^n for n in 4:10]
N_elements = [2^n for n in 2:8]

# allocation
l2e = zeros(length(N_datapoints),length(N_elements))

for (i,N_d) in enumerate(N_datapoints), (j,N_e) in enumerate(N_elements)
    @show N_d, N_e
    local dataset = create_dataset(N_d, x -> bar_E * x, strain_limit[2], strain_limit[1])

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
                                    force_func = force_func                                    
                                )

    local result = Datasolver.directSolverNonLinearBarA(
                        initProblem=initProblem,
                        constrained_dofs_global=constrained_dofs,
                        externalForce = (x,λ) -> force_func(α, β, x, λ),
                        dataset=dataset,
                        num_load_steps=num_load_steps,
                        loadFac=Vector(loadFac),
                        verbose=true,
                        NR_max_iter=20
                    );

    l2e[i,j] = Datasolver.relL2err1D(problem=initProblem, uNodal=result.u[end], uAfunction=x->uRef(x,β))
end

# plots

lvs = collect(-4:0.15:-2)
lvs_labels = map(x -> "$x", lvs)

contour(
        N_elements, N_datapoints, log10.(l2e),
        ylabel = "Number of data points",
        xlabel = "Number of elements",
        #levels = lvs,
        #levels_label=lvs_labels,
        clabels=false, 
        color=:jet,
        fill=false,
        colorbar=true,
        linewidth=2,
        #colorbar_title = L"Relative $L^2$ error (log10)",
        #colorbar_orientation = :horizontal,
        scale = :log10,
        framestyle = :box,
        dpi=150, 
        size=(1600,1200), 
        tickfont=font(16), 
        guidefont=font(16),
        legendfont=font(18)
)


savefig("fig/1DnonlinBar_convergenceADM.png")
savefig("fig/1DlinBar_convergenceADM.png")



# GoADM
l2eGA = zeros(length(N_datapoints),length(N_elements))

for (i,N_d) in enumerate(N_datapoints), (j,N_e) in enumerate(N_elements)
#for i in 3:8, j = [5,6,7]

    N_d = N_datapoints[i]
    N_e = N_elements[j]

    @show N_d, N_e
    
    local dataset = create_dataset(N_d, x -> bar_E * x, strain_limit[2], strain_limit[1])

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
                        num_load_steps=num_load_steps,
                        loadFac=Vector(loadFac),
                        verbose=true,
                        NR_max_iter=20
                    );

    l2eGA[i,j] = Datasolver.relL2err1D(problem=initProblem, uNodal=result2.u[end], uAfunction=x->uRef(x,β))
end



lvs = collect(-4:0.15:-2)
lvs_labels = map(x -> "$x", lvs)

contour(
        N_elements, N_datapoints, log10.(l2eGA),
        ylabel = "Number of data points",
        xlabel = "Number of elements",
        #levels = lvs,
        #levels_label=lvs_labels,
        clabels=false, 
        color=:jet,
        fill=false,
        colorbar=true,
        linewidth=2,
        #colorbar_title = L"Relative $L^2$ error (log10)",
        #colorbar_orientation = :horizontal,
        scale = :log10,
        framestyle = :box,
        dpi=150, 
        size=(1600,1200), 
        tickfont=font(16), 
        guidefont=font(16),
        legendfont=font(18) 
)



savefig("fig/1DnonlinBar_convergenceGoADM.png")
savefig("fig/1DlinBar_convergenceGoADM.png")



