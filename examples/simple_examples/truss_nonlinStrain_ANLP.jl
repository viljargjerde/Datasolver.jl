
using Datasolver, Revise, LinearAlgebra, Test, Plots, LaTeXStrings

include("examples/simple_examples/addFuncs.jl")

##### nonlinear (strain) truss structure - solving with mixed 3-field formulation ######


A = 2000/1e6        # [m²]
bar_E = sqrt(2)*1e+04   # [Pa] 

alpha = 1.0         # 0: linear strain     1: nonlinear strain

a = 1.0             # width and height of each bar

F = [0.0;4]         # [N]

# reference linear solution computed with https://valdivia.staff.jade-hs.de/fachwerk_en.html
uRef = [0;0;0;sqrt(2)*F[2]*a/A/bar_E;0;0]
eRef = [F[2]/A/bar_E/sqrt(2);F[2]/A/bar_E/sqrt(2)]
sRef = eRef .* bar_E

node_vector = [
    [0,   0],
    [a,   a],
    [2a,  0]
]

constrained_dofs = [
    (1, 1),
    (1, 2),
    (3, 1),
    (3, 2)
]

connections = [
    (1, 2),
    (2, 3)
]

force = zeros(2 * length(node_vector))
force[3:4] = F

initProblem = TrussProblem(
    A,
    force,
    connections,
    alpha,
    constrained_dofs,
    node_vector = node_vector,
    num_quad_pts = 2,
)

num_load_steps = 100

loadFac = LinRange(0.0,num_load_steps,num_load_steps+1)

nonlin_results = solveANLP(
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    YoungModulus = bar_E,
    NR_max_iter=20,
    qrFactorized = true
);


nonlin_results.NRiter

nonlin_results.u[end]

nonlin_results.e[end]

nonlin_results.s[end]

#--------------------------------------------------------

