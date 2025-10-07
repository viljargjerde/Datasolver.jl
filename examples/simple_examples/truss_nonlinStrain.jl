
using Datasolver, Revise, LinearAlgebra, Test, Plots, LaTeXStrings


##### nonlinear (strain) truss structure ######


A = 2000/1e6        # [m²]
bar_E = sqrt(2)*1e+04   # [Pa] 

alpha = 1.0         # 0: linear strain     1: nonlinear strain

NR_max_iter = 50
NR_tol = 1e-10

a = 1.0             # width and height of each bar

F = [0.0;4]         # [N]

numDataPts = 512

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


# strain limit = max of eRef + safety increase
strain_limit = 1.1 .* [eRef[1];
                       -eRef[1]]

dataset = create_dataset(numDataPts, x -> bar_E * x, strain_limit[2], strain_limit[1])


#-----------------------1 load step---------------------------------

    problem = TrussProblem(
        A,
        force,
        connections,
        alpha,
        constrained_dofs,
        node_vector = node_vector,
        num_quad_pts = 2,
    )

    nonlin_result = Datasolver.directSolverNonLinearBar(problem, dataset, NR_tol = NR_tol);

    uh = nonlin_result.u[end]
    eh = nonlin_result.e[end]
    sh = nonlin_result.s[end]

    uxNodal = uh[1:2:end]
    uyNodal = uh[2:2:end]

#--------------------------------------------------------




#-----------------------multiple load steps---------------------------------

initProblem = TrussProblem(
    A,
    force,
    connections,
    alpha,
    constrained_dofs,
    node_vector = node_vector,
    num_quad_pts = 2,
)

num_load_steps = 10

loadFac = LinRange(0,1.0,num_load_steps+1)

nonlin_results = Datasolver.directSolverNonLinearBarA(
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    dataset=dataset,
    verbose=true
);

nonlin_results.ADMiter'
nonlin_results.NRiter

sum(nonlin_results.ADMiter)

nonlin_results.u[end]

#--------------------------------------------------------

