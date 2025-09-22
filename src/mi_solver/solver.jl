export directSolverLinearBar, directSolverNonLinearBar
using Datasolver: SolveResults

function directSolverNonLinearBar(;
	node_vector::AbstractArray,
	data_set::AbstractArray,
	costFunc_ele::Function,
	random_init_data::Bool = true,
	DD_max_iter::Int = 2000,
	DD_tol::Float64 = 1e-10,
	num_ele::Int = 2,
	numQuadPts::Int = 2,
	costFunc_constant::Float64 = 1.0,
	bar_distF::Float64 = 1.0,
	cross_section_area::Float64 = 1.0,
	NR_num_load_step::Int = 50,
	NR_tol::Float64 = 1e-10,
	NR_max_iter::Int = 50,
)

	## initialize e_star and s_star
	numDataPts = size(data_set, 1)
	results = SolveResults(N_datapoints = numDataPts, Φ = node_vector)

	if random_init_data
		init_data_id = rand(1:numDataPts, num_ele)
		init_data = data_set[init_data_id, :]
	end

	data_star = deepcopy(init_data)

	## Newton_Raphson scheme combined in the direct solver [Kirchdoerfer - 2016 - Data-driven computational mechanics]
	# ndofs
	num_node = length(node_vector)
	ndof_u = ndof_lambda = num_node
	ndof_e = ndof_s = ndof_mu = num_ele

	ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
	ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]

	# boundary conditions: fixed-free
	constrained_dofs = [1 (ndof_u + ndof_e + ndof_s + ndof_mu + 1)]

	# allocation solution vector
	x = spzeros(ndof_tot)

	# iterative data-driven direct solver
	dd_iter = 0

	while dd_iter <= DD_max_iter

		# newton-raphson scheme
		x = spzeros(ndof_tot)        # initial guess for NR scheme

		for cc_load ∈ 1:NR_num_load_step
			iter = 0
			alpha = cc_load / NR_num_load_step

			while iter <= NR_max_iter
				iter += 1

				Delta_x = NewtonRaphsonStep(
					previous_sol = x,
					data_star = data_star,
					node_vector = node_vector,
					num_ele = num_ele,
					numQuadPts = numQuadPts,
					ndofs = ndofs,
					costFunc_constant = costFunc_constant,
					bar_distF = bar_distF * alpha,
					cross_section_area = cross_section_area,
					constrained_dofs = constrained_dofs,
				)

				# recover full dimension
				Delta_x = [0; Delta_x]
				Delta_x = [Delta_x[begin:ndof_u+ndof_e+ndof_s+ndof_mu]; 0; Delta_x[ndof_u+ndof_e+ndof_s+ndof_mu+1:end]]

				# update solution
				x += Delta_x

				# check convergence
				if norm(Delta_x) <= NR_tol
					break
				end


			end
			# @show iter

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
		data_star_new = assignLocalState(data_set = data_set, local_state = [ebar sbar], costFunc_ele = costFunc_ele)


		# evaluate global cost function (discrete)
		# push!(costFunc_global, integrateCostfunction(costFunc_ele = costFunc_ele, local_state = [ebar sbar], data_star = data_star, node_vector = node_vector, num_ele = num_ele, numQuadPts = numQuadPts, cross_section_area = cross_section_area))
		curr_cost = integrateCostfunction(costFunc_ele = costFunc_ele, local_state = [ebar sbar], data_star = data_star, node_vector = node_vector, num_ele = num_ele, numQuadPts = numQuadPts, cross_section_area = cross_section_area)

		## test convergence
		data_diff = data_star - data_star_new
		err = [norm(data_diff[i, :]) for i in 1:num_ele]
		max_err = maximum(err)

		dd_iter += 1

		if max_err <= DD_tol
			@show dd_iter
			break
		end


		# overwrite local state
		data_star = deepcopy(data_star_new)
		E = data_star[:, 1]
		S = data_star[:, 2]
		equilibrium = equilibrium_eq(uhat, bar_distF, node_vector, cross_section_area, sbar, num_ele)
		push!(results.u, collect(uhat))
		push!(results.e, collect(ebar))
		push!(results.s, collect(sbar))
		push!(results.λ, collect(λ))
		push!(results.μ, collect(μ))
		push!(results.E, collect(E))
		push!(results.S, collect(S))
		push!(results.cost, curr_cost)
		push!(results.balance, equilibrium)
		# TODO compatibility
	end
	@show results.balance
	return results
end

function equilibrium_eq(uhat, f, node_vector, cross_section_area, s, num_ele, numQuadPts = 2)
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)
	N_matrix, dN_matrix = constructBasisFunctionMatrixLinearLagrange(evalPts = quad_pts)
	R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts = quad_pts)
	equilibrium = zeros(num_ele + 1)

	for cc_ele ∈ 1:num_ele      # loop over elements    
		active_dofs_u = collect(cc_ele:cc_ele+1)
		active_dofs_s = cc_ele

		# jacobian for the integration
		xi0, xi1 = node_vector[cc_ele:cc_ele+1]
		J4int = (xi1 - xi0) / 2

		# jacobian for derivative
		J4deriv = (xi1 - xi0) / 2
		duh = (dN_matrix ./ J4deriv)' * uhat[active_dofs_u]
		@show active_dofs_s
		@show ((dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int) * (R_matrix * s[active_dofs_s]) * (1.0 .+ duh) - N_matrix * (quad_weights .* J4int) * f)
		equilibrium[active_dofs_u] += ((dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int) * (R_matrix * s[active_dofs_s]) * (1.0 .+ duh) - N_matrix * (quad_weights .* J4int) * f)

	end
	@show sum(equilibrium)
	@show equilibrium
	equilibrium

end

function directSolverLinearBar(;
	node_vector::AbstractArray,
	data_set::AbstractArray,
	costFunc_ele::Function,
	random_init_data::Bool = true,
	max_iter::Int = 2000,
	tol::Float64 = 1e-10,
	num_ele::Int = 2,
	numQuadPts::Int = 2,
	costFunc_constant::Float64 = 1.0,
	bar_distF::Float64 = 1.0,
	cross_section_area::Float64 = 1.0,
)

	## initialize e_star and s_star
	numDataPts = size(data_set, 1)

	if random_init_data
		init_data_id = rand(1:numDataPts, num_ele)
		init_data = data_set[init_data_id, :]
	end

	data_star = deepcopy(init_data)

	## direct solver [Kirchdoerfer - 2016 - Data-driven computational mechanics]
	# ndofs
	num_node = length(node_vector)
	results = SolveResults(N_datapoints = size(data_set, 1), Φ = node_vector)

	ndof_u = ndof_lambda = num_node
	ndof_e = ndof_s = ndof_mu = num_ele

	ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
	ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]

	# assembly system matrix
	A = assembleLinearSystemMatrix(node_vector = node_vector, num_ele = num_ele, ndofs = ndofs, costFunc_constant = costFunc_constant, cross_section_area = cross_section_area)

	# boundary conditions: fixed-free
	constrained_dofs = [1 (ndof_u + ndof_e + ndof_s + ndof_mu + 1)]
	ids = collect(1:size(A, 1))
	deleteat!(ids, constrained_dofs)
	A = A[ids, ids]

	# check condition number
	κ = cond(Matrix(A))
	@show κ

	# allocation solution vector
	x = spzeros(ndof_tot - length(constrained_dofs))
	costFunc_global = []

	# iterative direct solver
	iter = 0

	while iter <= max_iter
		# assembly rhs
		rhs = assembleRhsLinearBar(data_star = data_star, node_vector = node_vector, num_ele = num_ele, ndofs = ndofs, costFunc_constant = costFunc_constant, bar_distF = bar_distF, cross_section_area = cross_section_area)
		rhs = rhs[ids]

		# solving
		x = qr(Matrix(A)) \ rhs

		# check residual (Ax-rhs)
		r = norm(rhs - A * x)
		@show r
		full_x = [0; x]
		full_x = [full_x[1:ndof_u+ndof_e+ndof_s+ndof_mu]; 0; full_x[ndof_u+ndof_e+ndof_s+ndof_mu+1:end]]

		# collect computed ebar and sbar
		indices = cumsum(ndofs)

		# Extract variables from x using computed indices
		uhat = full_x[1:indices[1]]
		ebar = full_x[indices[1]+1:indices[2]]
		sbar = full_x[indices[2]+1:indices[3]]
		μ = full_x[indices[3]+1:indices[4]]
		λ = full_x[indices[4]+1:end]

		## local state assignment
		data_star_new = assignLocalState(data_set = data_set, local_state = [ebar sbar], costFunc_ele = costFunc_ele)

		# evaluate global cost function (discrete)
		curr_cost = integrateCostfunction(costFunc_ele = costFunc_ele, local_state = [ebar sbar], data_star = data_star, node_vector = node_vector, num_ele = num_ele, numQuadPts = numQuadPts, cross_section_area = cross_section_area)
		E = data_star_new[:, 1]
		S = data_star_new[:, 2]
		equilibrium = equilibrium_eq(uhat, bar_distF, node_vector, cross_section_area, sbar, num_ele)

		push!(results.u, collect(uhat))
		push!(results.e, collect(ebar))
		push!(results.s, collect(sbar))
		push!(results.λ, collect(λ))
		push!(results.μ, collect(μ))
		push!(results.E, collect(E))
		push!(results.S, collect(S))
		push!(results.cost, curr_cost)
		push!(results.balance, equilibrium)
		## test convergence
		data_diff = data_star - data_star_new
		err = [norm(data_diff[i, :]) for i in 1:num_ele]
		max_err = maximum(err)

		iter += 1

		if max_err <= tol
			@show iter
			break
		end

		# overwrite local state
		data_star = deepcopy(data_star_new)
	end

	# collect solution fields

	return results
end



function assignLocalState(; data_set::AbstractArray, local_state::AbstractArray, costFunc_ele::Function)
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



function integrateCostfunction(; costFunc_ele::Function, local_state::AbstractArray, data_star::AbstractArray, node_vector::AbstractArray, num_ele::Int = 2, numQuadPts::Int = 2, cross_section_area::Float64 = 1.0)

	# quad points in default interval [-1,1]
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)

	# integration
	costFunc_global = 0

	for e in 1:num_ele      # loop over element
		# jacobian for the integration
		xi0, xi1 = node_vector[e:e+1]
		J4int = (xi1 - xi0) / 2

		costFunc_global += costFunc_ele(local_state[e, 1] - data_star[e, 1], local_state[e, 2] - data_star[e, 2]) * sum(quad_weights) * J4int * cross_section_area
	end

	return costFunc_global
end
