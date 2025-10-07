
using Datasolver, Revise, LinearAlgebra, Test, Plots, LaTeXStrings


##### nonlinear (strain) 1D bar structure ######

bar_L = Float64(π)
A = 2000/1e6        # [m²]
bar_E = 1.622e+03   # [Pa] 

α = 1.0         # 0: linear strain     1: nonlinear strain

ne = 16

numDataPts = 1024

# manifactured solution and nonlinear force function
β = 0.1*π

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
strain_limit = 1.1 .* [maximum(x[1] for x in eRef.(xx,α,β));
                        minimum(x[1] for x in eRef.(xx,α,β))]

dataset = create_dataset(numDataPts, x -> bar_E * x, strain_limit[2], strain_limit[1])


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
        num_load_steps=num_load_steps,
        loadFac=Vector(loadFac),
        verbose=true
    );
  


# taking results of the last ADM iter
uh = nonlin_result.u[end]
uxh = uh[1:2:end]

eh = nonlin_result.e[end]
sh = nonlin_result.s[end]


# plots
plot(xx, uRef.(xx,β), linewidth=2, linecolor=:black,label="uRef")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], uxh, linewidth=2, linecolor=:royalblue,label="uh")


plot(xx,eRef.(xx,α,β), linewidth=2, linecolor=:black,label="eRef")
plot!([node_vector[i][1] for i in 1:initProblem.num_node], [eh[1];eh], linewidth=2,linetype=:steppre,label="eh")


plot(xx,sRef.(xx,α,β), linewidth=2, linecolor=:black)
plot!([node_vector[i][1] for i in 1:initProblem.num_node], [sh[1];sh], linewidth=2,linetype=:steppre)


plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="Displacement", tickfont=font(16), guidefont=font(16),legendfont=font(18))


savefig("fig/1DnonlinBar_uh.png")
savefig("fig/1DnonlinBar_eh.png")



#----------------------------------------------------------------
# CONVERGENCE STUDY
bar_L = Float64(π)
A = 2000/1e6        # [m²]
bar_E = 1.622e+03   # [Pa] 

α = 1.0
# manifactured solution and nonlinear force function
β = 0.1*π

uRef(x,β) = β * sin(π*x/bar_L)

eRef(x,α,β) = β * π / bar_L * cos(π*x/bar_L) * ( 1 + 0.5*α * β * π / bar_L * cos(π*x/bar_L) )

sRef(x,α,β) = bar_E * eRef(x,α,β)

force_func(α, β, x, λ) = λ .*
    [bar_E * A * ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_L))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_L)) / bar_L))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_L)) / bar_L) + 2) / bar_L^(2)); 0]

num_load_steps = 5

loadFac = LinRange(0,1.0,num_load_steps+1)

# loop over num_ele and numDataPts

N_datapoints = [2^n for n in 2:9]
N_elements = [2^n for n in 2:9]

xx = 0:bar_L/1000:bar_L
strain_limit = 1.1 .* [maximum(x[1] for x in eRef.(xx,α,β));
                    minimum(x[1] for x in eRef.(xx,α,β))]

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
                        verbose=true
                    );

    l2e[i,j] = Datasolver.relL2err1D(problem=initProblem, uNodal=result.u[end], uAfunction=x->uRef(x,β))
end

# plots
contour(
        N_elements, N_datapoints, log10.(l2e'),
        ylabel = "Number of data points",
        xlabel = "Number of elements",
        colorbar_title = L"Relative $L^2$ error (log10)",
        scale = :log10,
        fill = false,
        framestyle = :box,
        dpi=150, 
        size=(800,600), 
        tickfont=font(16), 
        guidefont=font(16),
        legendfont=font(18) )

savefig("fig/1DnonlinBar_convergence.png")

############## linear strain

α = 0.0

strain_limit = 1.1 .* [maximum(x[1] for x in eRef.(xx,α,β));
                    minimum(x[1] for x in eRef.(xx,α,β))]

# allocation
l2e_lin = zeros(length(N_datapoints),length(N_elements))

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
                        verbose=true
                    );

    l2e_lin[i,j] = Datasolver.relL2err1D(problem=initProblem, uNodal=result.u[end], uAfunction=x->uRef(x,β))
end

contour(
        N_elements, N_datapoints, log10.(l2e_lin'),
        ylabel = "Number of data points",
        xlabel = "Number of elements",
        colorbar_title = L"Relative $L^2$ error (log10)",
        scale = :log10,
        fill = false,
        framestyle = :box,
        dpi=150, 
        size=(800,600), 
        tickfont=font(16), 
        guidefont=font(16),
        legendfont=font(18) )

savefig("fig/1DlinBar_convergence.png")
