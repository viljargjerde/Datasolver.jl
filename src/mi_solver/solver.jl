export directSolverLinearBar, directSolverNonLinearBar
using Datasolver: SolveResults
using NonlinearSolve
function directSolverNonLinearBar(
	node_vector::AbstractArray,
	constrained_dofs::AbstractArray,
	data_set::AbstractArray,
	costFunc_ele::Function,
	cross_section_area::Float64,
	bar_distF::Vector{Float64},
	num_ele::Int,
	costFunc_constant::Float64;
	random_init_data::Bool = true,
	DD_max_iter::Int = 2000,
	DD_tol::Float64 = 1e-10,
	numQuadPts::Int = 2,
	NR_num_load_step::Int = 50,
	NR_tol::Float64 = 1e-10,
	NR_max_iter::Int = 50,
	alpha::Float64 = 1.0,
	NR_damping::Float64 = 1.0,
)

	## initialize e_star and s_star
	numDataPts = size(data_set, 1)
	results = SolveResults(N_datapoints = numDataPts, Φ = node_vector)
	num_node = length(node_vector)
	dims = length(node_vector[1])
	ndof_u = ndof_lambda = num_node * dims
	ndof_e = ndof_s = ndof_mu = num_ele

	ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
	ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]

	if random_init_data
		init_data_id = rand(1:numDataPts, num_ele)
		init_data = data_set[init_data_id, :]
	else
		init_data = zeros(num_ele, 2)
		s = get_initialization_s(bar_distF, node_vector, cross_section_area, num_ele, ndofs, constrained_dofs, numQuadPts)
		S = zeros(length(s))
		E = zeros(length(s))
		for i in eachindex(S)
			# choose E, S pair where S is closest to s
			best_idx = argmin((abs(data_set[j, 2] - s[i]) for j in eachindex(data_set[:, 1])))
			init_data[i, :] = data_set[best_idx, :]
		end
	end

	data_star = deepcopy(init_data)

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
					data_star,
					node_vector,
					num_ele,
					ndofs,
					costFunc_constant,
					bar_distF * load_alpha,
					cross_section_area,
					alpha,
					constrained_dofs;
					numQuadPts = numQuadPts,
				)

				# recover full dimension. This can be optimized by removing dofs here instead of in the NR step, removing the need to reconstruct_vector in loop
				Delta_x = reconstruct_vector(Delta_x, constrained_dofs)

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
		data_star_new = assignLocalState(data_set, [ebar sbar], costFunc_ele)


		curr_cost = integrateCostfunction(costFunc_ele = costFunc_ele, local_state = [ebar sbar], data_star = data_star, node_vector = node_vector, num_ele = num_ele, numQuadPts = numQuadPts, cross_section_area = cross_section_area)

		## test convergence
		data_diff = data_star - data_star_new
		err = [norm(data_diff[i, :]) for i in 1:num_ele]
		max_err = maximum(err)

		dd_iter += 1

		# overwrite local state
		data_star = deepcopy(data_star_new)
		E = data_star[:, 1]
		S = data_star[:, 2]
		equilibrium = equilibrium_eq(uhat, bar_distF, node_vector, cross_section_area, sbar, num_ele, alpha, constrained_dofs)
		compat = compatibility_eq(uhat, node_vector, cross_section_area, ebar, num_ele, alpha)
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

		if max_err <= DD_tol
			break
		end
	end
	return results
end

function equilibrium_eq(uhat, f, node_vector, cross_section_area, s, num_ele, alpha, constrained_dofs, numQuadPts = 2)
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)
	dims = size(node_vector[1])[1]
	N_func, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)
	# R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts = quad_pts)
	equilibrium = zeros((num_ele + 1) * dims)

	for cc_ele ∈ 1:num_ele      # loop over elements    
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)

			active_dofs_u = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)

			active_dofs_s = cc_ele
			sh = s[active_dofs_s]

			# jacobian for the integration
			xi0, xi1 = node_vector[cc_ele:cc_ele+1]
			J4int = norm(xi1 - xi0) / 2

			# jacobian for derivative
			J4deriv = norm(xi1 - xi0) / 2
			dN_matrix = dN_func(quad_pt) / J4deriv
			N_matrix = N_func(quad_pt)
			duh = dN_matrix * uhat[active_dofs_u]
			dPhih = dN_matrix * [xi0; xi1]

			PBh = (dPhih + alpha * duh)
			integration_factor = cross_section_area * quad_weight * J4int

			equilibrium[active_dofs_u] += N_matrix' * (quad_weight * J4int * f) -
										  (dN_matrix') * (integration_factor) * (PBh * sh)

		end
	end
	idxs = collect(1:length(equilibrium))
	deleteat!(idxs, constrained_dofs[begin:length(constrained_dofs)÷2])
	equilibrium[idxs] # Remove constrained dofs

end


function compatibility_eq(uhat, node_vector, cross_section_area, e, num_ele, alpha, numQuadPts = 2)
	dims = size(node_vector[1])[1]
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)
	_, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)
	compatibility = zeros(num_ele)

	for cc_ele ∈ 1:num_ele      # loop over elements    
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)

			xi0, xi1 = node_vector[cc_ele:cc_ele+1]
			# jacobian for the integration
			J4int = norm(xi1 - xi0) / 2
			# jacobian for derivative
			J4deriv = norm(xi1 - xi0) / 2
			dN_matrix = dN_func(quad_pt) / J4deriv

			active_dofs_u = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
			active_dofs_e = cc_ele
			eh = e[active_dofs_e]
			duh = dN_matrix * uhat[active_dofs_u]
			integration_factor = cross_section_area * quad_weight * J4int
			dPhih = dN_matrix * [xi0; xi1]

			e_uh = duh' * dPhih + alpha / 2 * duh' * duh

			compatibility[active_dofs_e] += -(integration_factor * (e_uh - eh))
		end
	end
	compatibility

end



function assignLocalState(data_set::AbstractArray, local_state::AbstractArray, costFunc_ele::Function)
	num_local_state = size(local_state, 1)
	numDataPts = size(data_set, 1)

	# allocation
	data_star = zeros(num_local_state, 2)

	# find the closest data point to the local state
	for ii ∈ 1:num_local_state
		d = zeros(numDataPts)

		for i ∈ 1:numDataPts
			d[i] = sqrt(costFunc_ele(local_state[ii, 1] - data_set[i, 1], local_state[ii, 2] - data_set[i, 2]))
		end
		min_id = findmin(d)[2]
		data_star[ii, :] = data_set[min_id, :]
	end

	return data_star
end



function get_initialization_s(f, node_vector, cross_section_area, num_ele, ndofs, constrained_dofs, numQuadPts = 2)
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)
	dims = size(node_vector[1], 1)
	N_func, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)
	# R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts = quad_pts)
	ndof_u, _, ndof_s, _, _ = ndofs

	# Construct global matrices A sbar = b
	A = spzeros(ndof_u, ndof_s)
	b = spzeros(ndof_u)
	for cc_ele ∈ 1:num_ele      # loop over elements    
		active_dofs_u = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
		active_dofs_s = cc_ele
		# jacobian for the integration
		xi0, xi1 = node_vector[cc_ele:cc_ele+1]
		J4int = norm(xi1 - xi0) / 2

		# jacobian for derivative
		J4deriv = norm(xi1 - xi0) / 2
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)
			dN_matrix = dN_func(quad_pt) / J4deriv
			N_matrix = N_func(quad_pt)
			dPhih = dN_matrix * [xi0; xi1]

			integration_factor = cross_section_area * quad_weight * J4int
			A[active_dofs_u, active_dofs_s] += integration_factor * dN_matrix' * dPhih
			b[active_dofs_u] += N_matrix' * (quad_weight * J4int * f)

		end

	end
	idxs = collect(1:length(b))
	deleteat!(idxs, constrained_dofs[begin:length(constrained_dofs)÷2]) # constrained_dofs include lambda, but here we only care about u, which is the first half
	Matrix(A[idxs, :]) \ b[idxs]
end


function integrateCostfunction(; costFunc_ele::Function, local_state::AbstractArray, data_star::AbstractArray, node_vector::AbstractArray, num_ele::Int, cross_section_area::Float64, numQuadPts::Int = 2)

	# quad points in default interval [-1,1]
	_, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)

	# integration
	costFunc_global = 0

	for e in 1:num_ele      # loop over element
		# jacobian for the integration
		xi0, xi1 = node_vector[e:e+1]
		J4int = norm(xi1 - xi0) / 2

		costFunc_global += costFunc_ele(local_state[e, 1] - data_star[e, 1], local_state[e, 2] - data_star[e, 2]) * sum(quad_weights) * J4int * cross_section_area
	end

	return costFunc_global
end
