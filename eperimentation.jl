##
using Revise
using Datasolver
using Symbolics
using LinearAlgebra
L = 1.0
E = 2e5
data = create_dataset(500, x -> E * x, -15.0, 15.0)
A = 0.002
@variables x
connections, Φ, f, _ = setup_1d_bar_problem(2, L, x -> x)
f = 1000.0
Φ = [p[1] for p in Φ]
fixed_dof = [1]
##
Np = Datasolver.construct_Np(Φ)
N = Datasolver.construct_N(Φ)
R = Datasolver.construct_R(Φ, connections)
m = length(Np) # hat / nodes 
n = length(R) # bar / elements
row_sizes = [m, n, n, n, m]
col_sizes = [m, n, n, n, m]
u = zeros(m)
s = zeros(n)
A_mat = Datasolver.A(u, Np, R, s, data.C, A, L, connections, Φ)
J_mat = Datasolver.J(zeros(sum(row_sizes)), Np, R, data.C, A, L, connections, Φ)
@show fixed_dof

##

Datasolver.nonlin_datasolve(connections, Φ, A, data, f, L, fixed_dof; initialization = true)

## 
extract_A(row, col) = Datasolver.extract_matrix_block(A_mat, row_sizes, col_sizes, row, col)
mat = [
	extract_A(1, 4) extract_A(1, 5)
	extract_A(2, 4) zeros(n, m)
	zeros(n, n) extract_A(3, 5)
]
@show rank(mat), minimum(size(mat))

## LP - Kanno

L = 3.6 # m
A = 0.002
F = 4000 # N
N_datapoints = 2000
strain_magnitude = 0.005
noise_magnitude = 0.0

connections = [[1, 3], [3, 5], [3, 4], [5, 6], [2, 4], [4, 6], [1, 4], [2, 3], [3, 6], [4, 5]]
Φ = [[0, 0], [0, L], [L, 0], [L, L], [2L, 0], [2L, L]]
dataset = create_dataset(N_datapoints, x -> 5e6 * tanh.(500 .* x), -strain_magnitude, strain_magnitude, noise_magnitude = noise_magnitude)
fixed_dof = [(1, 1), (1, 2), (2, 1), (2, 2)]
f = zeros(2 * length(Φ) - 4)
f[2] = -F
f[6] = -F
result = Datasolver.LP_solver(connections, Φ, A, dataset, f, fixed_dof, verbose = false)
plot_dataset(dataset, Datasolver.get_final(result))



## LP 1d bar 

L = 1.0
A = 0.002
E = 2e5
N_elements = 20
N_datapoints = 100
E_3 = 5000.0
dataset = create_dataset(N_datapoints, x -> E * x + E_3 * x^3, -5.0, 5.0; noise_magnitude = 0.0)
f_0 = 4000.0
n = 1
connections, Φ, f, fixed_dof = setup_1d_bar_problem(N_elements, L, x -> f_0 * sin(n * pi * x / L))
result = Datasolver.LP_solver(connections, Φ, A, dataset, f, fixed_dof; verbose = false)

# plot_dataset(dataset, Datasolver.get_final(result);title=L"σ = E_1 ϵ + E_3 ϵ^3")
plot_results(result, dataset = dataset)

##

