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
A_mat = Datasolver.A(u, Np, R, s, data.C, A, connections, Φ)
J_mat = Datasolver.J(zeros(sum(row_sizes)), Np, R, data.C, A, connections, Φ)
@show fixed_dof

##

Datasolver.nonlin_datasolve(connections, Φ, A, data, f, L, fixed_dof; initialization = true)

## 

### This is still not full rank ###

extract_A(row, col) = Datasolver.extract_matrix_block(A_mat, row_sizes, col_sizes, row, col)
mat = [
	extract_A(1, 4) extract_A(1, 5)
	extract_A(2, 4) zeros(n, m)
	zeros(n, n) extract_A(3, 5)
]
@show rank(mat), minimum(size(mat))
## 
