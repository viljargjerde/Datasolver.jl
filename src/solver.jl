using LinearAlgebra


costFunc_ele = (e, s, C) -> 0.5 * (C * e^2 + 1 / C * s^2);


"""
	check_dataset_is_safe(E, S, data) -> Bool

Checks if `E` and `S` are within safe boundaries (not touching the dataset edges).

# Arguments
- `E::Vector`: Strain values to check.
- `S::Vector`: Stress values to check.
- `data::Dataset`: Dataset containing strain and stress limits.

# Returns
- `true` if the values are within safe boundaries, otherwise `false`.
"""
function check_dataset_is_safe(E, S, data)
	for e in E
		if e == data.E[begin] || e == data.E[end]
			return false
		end
	end
	for s in S
		if s == data.S[begin] || s == data.S[end]
			return false
		end
	end
	return true
end


function directSolverNonLinearBar(
	problem::Barproblem,
	dataset::Dataset;
	random_init_data::Bool = false,
	DD_max_iter::Int = 100,
	NR_num_load_step::Int = 50,
	NR_tol::Float64 = 1e-10,
	NR_max_iter::Int = 50,
	NR_damping::Float64 = 1.0,
)

	## initialize e_star and s_star
	numDataPts = length(dataset)
	node_vector = problem.node_vector
	num_ele = problem.num_ele
	results = SolveResults(N_datapoints = numDataPts, Φ = node_vector)
	num_node = length(node_vector)
	dims = length(node_vector[1])
	ndof_u = ndof_lambda = num_node * dims
	ndof_e = ndof_s = ndof_mu = num_ele

	ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
	ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]

	if random_init_data
		init_data_id = rand(1:numDataPts, num_ele)
		E, S = dataset[init_data_id]
	else
		s = get_initialization_s(problem)
		S = zeros(length(s))
		E = zeros(length(s))
		for i in eachindex(S)
			# choose E, S pair where S is closest to s
			best_idx = argmin((abs(dataset.S[j] - s[i]) for j in eachindex(dataset.S)))
			S[i] = dataset.S[best_idx]
			E[i] = dataset.E[best_idx]
		end
	end

	x = spzeros(ndof_tot)

	# iterative data-driven direct solver
	dd_iter = 0

	while dd_iter <= DD_max_iter

		# newton-raphson scheme
		# x = spzeros(ndof_tot)
		x[ndof_u+ndof_e+ndof_s+1:ndof_u+ndof_e+ndof_s+ndof_mu] = rand(ndof_mu)
		for cc_load ∈ 1:NR_num_load_step
			iter = 0
			load_alpha = cc_load / NR_num_load_step
			# x[ndof_u+ndof_e+ndof_s+1:ndof_u+ndof_e+ndof_s+ndof_mu] = rand(ndof_mu)

			while iter <= NR_max_iter
				iter += 1

				Delta_x = NewtonRaphsonStep(
					x,
					E,
					S,
					dataset.C,
					load_alpha,
					problem,
				)

				# recover full dimension. This can be optimized by removing dofs here instead of in the NR step, removing the need to reconstruct_vector in loop
				Delta_x = reconstruct_vector(Delta_x, problem.constrained_dofs)

				# update solution
				x += NR_damping * Delta_x

				# check convergence
				if norm(Delta_x) <= NR_tol
					break
				end


			end

			if iter == NR_max_iter
				println("NR did not converge at load step ")
				@show cc_load
				break
			end
		end



		# collect computed ebar and sbar
		indices = cumsum(ndofs)

		# Extract variables from x using computed indices
		uhat = x[1:indices[1]]
		ebar = x[indices[1]+1:indices[2]]
		sbar = x[indices[2]+1:indices[3]]
		μ = x[indices[3]+1:indices[4]]
		λ = x[indices[4]+1:end]
		## local state assignment
		(new_E, new_S) = assignLocalState(dataset, ebar, sbar)


		curr_cost = integrateCostfunction(ebar, sbar, E, S, dataset.C, problem)

		converged = (new_E == E) && (new_S == S)
		dd_iter += 1

		# overwrite local state
		E = new_E
		S = new_S
		equilibrium = equilibrium_eq(uhat, sbar, problem)
		compat = compatibility_eq(uhat, ebar, problem)
		push!(results.u, [norm(uhat[i:i+dims-1]) for i in 1:dims:length(uhat)])
		push!(results.e, collect(ebar))
		push!(results.s, collect(sbar))
		push!(results.λ, [norm(λ[i:i+dims-1]) for i in 1:dims:length(λ)])
		push!(results.μ, collect(μ))
		push!(results.E, collect(E))
		push!(results.S, collect(S))
		push!(results.cost, curr_cost)
		push!(results.equilibrium, equilibrium)
		push!(results.compatibility, compat)

		if converged
			break
		end
	end
	return results
end

function equilibrium_eq(uhat, sbar, problem::Barproblem)
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)
	dims = problem.dims
	N_func, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)
	equilibrium = zeros((problem.num_node) * dims)
	alpha = problem.alpha
	for cc_ele ∈ 1:problem.num_ele      # loop over elements    
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)

			active_dofs_u = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)

			active_dofs_s = cc_ele
			sh = sbar[active_dofs_s]

			# jacobian for the integration
			xi0, xi1 = problem.node_vector[cc_ele:cc_ele+1]
			J4int = norm(xi1 - xi0) / 2

			# jacobian for derivative
			J4deriv = norm(xi1 - xi0) / 2
			dN_matrix = dN_func(quad_pt) / J4deriv
			N_matrix = N_func(quad_pt)
			duh = dN_matrix * uhat[active_dofs_u]
			dPhih = dN_matrix * [xi0; xi1]

			PBh = (dPhih + alpha * duh)
			integration_factor = problem.area * quad_weight * J4int

			equilibrium[active_dofs_u] += N_matrix' * (quad_weight * J4int * problem.force(quad_pt)) -
										  (dN_matrix') * (integration_factor) * (PBh * sh)

		end
	end
	idxs = collect(1:length(equilibrium))
	deleteat!(idxs, problem.constrained_dofs[begin:length(problem.constrained_dofs)÷2])
	equilibrium[idxs] # Remove constrained dofs

end


function compatibility_eq(uhat, ebar, problem::Barproblem)
	dims = problem.dims
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)
	_, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)
	compatibility = zeros(problem.num_ele)
	alpha = problem.alpha

	for cc_ele ∈ 1:problem.num_ele      # loop over elements    
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)

			xi0, xi1 = problem.node_vector[cc_ele:cc_ele+1]
			# jacobian for the integration
			J4int = norm(xi1 - xi0) / 2
			# jacobian for derivative
			J4deriv = norm(xi1 - xi0) / 2
			dN_matrix = dN_func(quad_pt) / J4deriv

			active_dofs_u = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
			active_dofs_e = cc_ele
			eh = ebar[active_dofs_e]
			duh = dN_matrix * uhat[active_dofs_u]
			integration_factor = problem.area * quad_weight * J4int
			dPhih = dN_matrix * [xi0; xi1]

			e_uh = duh' * dPhih + alpha / 2 * duh' * duh

			compatibility[active_dofs_e] += -(integration_factor * (e_uh - eh))
		end
	end
	compatibility

end



function assignLocalState(dataset::Dataset, ebar::AbstractArray, sbar::AbstractArray)

	# allocation
	E = zero(ebar)
	S = zero(sbar)

	# find the closest data point to the local state
	for i ∈ eachindex(E)

		distances = (sqrt(costFunc_ele(dataset.E[j] - ebar[i], dataset.S[j] - sbar[i], dataset.C)) for j in 1:length(dataset))
		min_id = argmin(distances)
		E[i], S[i] = dataset[min_id]
	end

	return (E, S)
end



function get_initialization_s(problem::Barproblem)
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)
	dims = problem.dims
	N_func, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)
	ndof_u, _, ndof_s, _, _ = get_ndofs(problem)

	# Construct global matrices A sbar = b
	A = spzeros(ndof_u, ndof_s)
	b = spzeros(ndof_u)
	for cc_ele ∈ 1:problem.num_ele      # loop over elements    
		active_dofs_u = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
		active_dofs_s = cc_ele
		# jacobian for the integration
		xi0, xi1 = problem.node_vector[cc_ele:cc_ele+1]
		J4int = norm(xi1 - xi0) / 2

		# jacobian for derivative
		J4deriv = norm(xi1 - xi0) / 2
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)
			dN_matrix = dN_func(quad_pt) / J4deriv
			N_matrix = N_func(quad_pt)
			dPhih = dN_matrix * [xi0; xi1]

			integration_factor = problem.area * quad_weight * J4int
			A[active_dofs_u, active_dofs_s] += integration_factor * dN_matrix' * dPhih
			b[active_dofs_u] += N_matrix' * (quad_weight * J4int * problem.force(quad_pt))

		end

	end
	idxs = collect(1:length(b))
	deleteat!(idxs, problem.constrained_dofs[begin:length(problem.constrained_dofs)÷2]) # constrained_dofs include lambda, but here we only care about u, which is the first half
	Matrix(A[idxs, :]) \ b[idxs]
end


function integrateCostfunction(e::AbstractArray, s::AbstractArray, E::AbstractArray, S::AbstractArray, costFunc_constant::Float64, problem::Barproblem)

	# quad points in default interval [-1,1]
	_, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)

	# integration
	costFunc_global = 0

	for i in 1:problem.num_ele      # loop over element
		# jacobian for the integration
		xi0, xi1 = problem.node_vector[i:i+1]
		J4int = norm(xi1 - xi0) / 2

		costFunc_global += costFunc_ele(e[i] - E[i], s[i] - S[i], costFunc_constant) * sum(quad_weights) * J4int * problem.area
	end

	return costFunc_global
end



function NewtonRaphsonStep(
	x::AbstractArray,
	E::AbstractArray,
	S::AbstractArray,
	costFunc_constant,
	load_alpha::Float64,
	problem::Barproblem,
)

	# assembly
	rhs = assembleEquilibriumResidual(
		x,
		E,
		S,
		costFunc_constant,
		problem,
		load_alpha,
	)

	J = assembleLinearizedSystemMatrix(x, problem, costFunc_constant)

	# enforcing boundary conditions    
	ids = collect(1:size(J, 1))
	deleteat!(ids, problem.constrained_dofs)
	J = J[ids, ids]
	rhs = rhs[ids]
	# solving
	Delta_x = qr(Matrix(J)) \ rhs

	# check residual and condition number of J
	r = norm(rhs - J * Delta_x)
	κ = cond(Matrix(J))

	if r > 1e-10
		println("Warning: Solution with residual $r > 1e-10")
	end
	if κ > 1e20
		# @show κ
		println("Warning: Condition number of the system matrix $κ > 1e20")
	end

	return Delta_x
end



