using Datasolver

node_vector = [
    [0,     0],
    [3.6,   0],
    [2*3.6, 0],
    [2*3.6, 3.6],
    [3.6,   3.6],
    [0,     3.6]
]  
A = 2000/1e6       # [m²]

constrained_dofs = [
    (1, 1),
    (1, 2),
    (6, 1),
    (6, 2)
]

alpha = 0.0


connections = [
    (1, 2),
    (2, 3),
    (1, 5),
    (2, 6),
    (2, 5),
    (2, 4),
    (3, 5),
    (3, 4),
    (5, 6),
    (4, 5)
]
num_data_pts = 20
bar_E = 1.7446e+09

force = zeros(2 * length(node_vector))
λ = 1

force[4] = -400.0*λ      # [N]   - downward force at node 2
force[6] = -400.0*λ   # [N]   - downward force at node 3



problem = TrussProblem(
    A,
    force,  # [N]   - downward force at node 3
    connections,
    alpha,
    constrained_dofs,
    node_vector = node_vector,
    num_quad_pts = 2,
)


dataset = create_dataset(num_data_pts, x -> bar_E * x, -5e-3, 5e-3)

results = directSolverNonLinearBar(
    problem,
    dataset;
    random_init_data = false,
    NR_max_iter = 1000,
    NR_tol = 1e-4,
);

