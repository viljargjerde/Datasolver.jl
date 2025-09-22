
function setup_comp_case()
	L = 7.0
	A = 0.013
	E = 3e5
	data = create_dataset(200, x -> E * x, -15.0, 15.0)
	f0 = 1700.0
	connections, Φ, f_vec, fixed_dof = setup_1d_bar_problem(11, L, x -> f0, false)
	Φ = [p[1] for p in Φ]
	Np = Datasolver.construct_Np(Φ)
	N = Datasolver.construct_N(Φ)
	R = Datasolver.construct_R(Φ, connections)
	B = Datasolver.create_B_matrix_1D(connections, Φ)

	return L, A, E, data, connections, Φ, f_vec, f0, fixed_dof, Np, N, R, B

end

# Test for the Dataset struct
@testset "All block entries same" begin
	L, A, E, data, connections, Φ, f_vec, f, fixed_dof, Np, N, R, B = setup_comp_case()
	m = length(Np) # hat / nodes 
	n = length(R) # bar / elements
	u = zeros(m)
	s = zeros(n)
	x_sizes = [m, n, n, n, m]
	x = vcat([zeros(s) for s in x_sizes]...)
	J_mat = Datasolver.J(x, Np, R, data.C, A, L, connections, Φ)
	w = Datasolver.calc_w(connections, Φ, A)
	C = data.C
	Z = zeros
	# TODO extract into function and test this function
	A = [
		Z(m, m) Z(m, n) Z(m, n) B' Z(m, m);
		Z(n, m) Diagonal(w .* C) Z(n, n) I Z(n, m);
		Z(n, m) Z(n, n) Diagonal(w ./ C) Z(n, n) B;
		-B I Z(n, n) Z(n, n) Z(n, m);
		Z(m, m) Z(m, n) (w .* B)' Z(m, n) Z(m, m)
	]

	extract(mat, row, col) = Datasolver.extract_matrix_block(mat, x_sizes, x_sizes, row, col)
	@testset "row: $r col $c" for r in 1:5, c in 1:5
		@test extract(A, r, c) ≈ extract(J_mat, r, c)

	end

end
@testset "Remove dofs" begin
	L, A, E, data, connections, Φ, f_vec, f, fixed_dof, Np, N, R, B = setup_comp_case()
	m = length(Np) # hat / nodes 
	n = length(R) # bar / elements
	u = zeros(m)
	s = zeros(n)
	x_sizes = [m, n, n, n, m]
	x = zeros(sum(x_sizes))
	J_mat = Datasolver.J(x, Np, R, data.C, A, L, connections, Φ)
	w = Datasolver.calc_w(connections, Φ, A)
	C = data.C
	Z = zeros

	n, m = size(B)

	A_mat = [
		Z(m, m) Z(m, n) Z(m, n) B' Z(m, m);
		Z(n, m) Diagonal(w .* C) Z(n, n) I Z(n, m);
		Z(n, m) Z(n, n) Diagonal(w ./ C) Z(n, n) B;
		-B I Z(n, n) Z(n, n) Z(n, m);
		Z(m, m) Z(m, n) (w .* B)' Z(m, n) Z(m, m)
	]

	@test J_mat ≈ A_mat

	B = B[:, 2:end]
	n, m = size(B)
	A_mat = [
		Z(m, m) Z(m, n) Z(m, n) B' Z(m, m);
		Z(n, m) Diagonal(w .* C) Z(n, n) I Z(n, m);
		Z(n, m) Z(n, n) Diagonal(w ./ C) Z(n, n) B;
		-B I Z(n, n) Z(n, n) Z(n, m);
		Z(m, m) Z(m, n) (w .* B)' Z(m, n) Z(m, m)
	]
	fixed_dofs = [1]
	free_dofs = Datasolver.get_free_dofs(fixed_dofs, x_sizes)
	J_free = J_mat[free_dofs, free_dofs]
	for i in eachindex(A_mat)
		if !(A_mat[i] ≈ J_free[i])
			@show i, A_mat[i] J_free[i]
		end
	end
	@test J_free ≈ A_mat

end

@testset "Δx is the same 1st iteration" begin
	L, A, E, data, connections, Φ, f_vec, f, fixed_dof, Np, N, R, B = setup_comp_case()
	B = B[:, 2:end]
	f_vec = f_vec[2:end]
	fixed_dof = [1]
	m = length(Np) # hat / nodes 
	n = length(R) # bar / elements
	x_sizes = [m, n, n, n, m]
	# Choose datapoints
	points = (data[rand(1:length(data))] for _ in 1:n)
	E, S = Datasolver.unpack_pairs(points)
	free_dofs = Datasolver.get_free_dofs(fixed_dof, x_sizes)

	w = Datasolver.calc_w(connections, Φ, A)
	e, s, u, λ, μ = Datasolver.solve_system(B, data.C, w, E, S, 2 * f_vec)
	x_lin = vcat(u, e, s, μ, λ)
	x0 = zeros(sum(x_sizes))
	J_mat = Datasolver.J(x0, Np, R, data.C, A, L, connections, Φ)
	g_vec = Datasolver.g(x0, N, Np, R, data.C, E, S, f, A, L, connections, Φ)
	Δx = -J_mat[free_dofs, free_dofs] \ g_vec[free_dofs]

	@test Δx ≈ x_lin
end

@testset "b is the same" begin
	L, A, E, data, connections, Φ, f_vec, f, fixed_dof, Np, N, R, B = setup_comp_case()
	f_vec = f_vec[2:end]
	w = Datasolver.calc_w(connections, Φ, A)
	Z = zeros
	int_L(integrand) = A * Datasolver.int_L_(connections, Φ, integrand)
	int_N(integrand) = 2 * Datasolver.int_N_(N, connections, Φ, integrand)
	k, n = size(B)
	C = data.C

	E = [data[i+length(data)÷2][1] for i in 1:k]
	S = [data[i+length(data)÷2][2] for i in 1:k]
	b_lin = vcat(
		Z(n),
		w .* C .* E,
		(w ./ C) .* S,
		Z(k),
		2 * f_vec, #! * 2??
	)



	m = length(Np)
	n = length(R)
	b_nl = [Z(m)
		int_L(C) * E
		int_L(1 / C) * S
		Z(n)
		int_N(f)[2:end]]
	@test b_lin ≈ b_nl


end;
