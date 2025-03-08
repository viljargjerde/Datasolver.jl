using LinearAlgebra

"""
	struct Dataset{T <: Real}

Represents a dataset containing vectors of empirical strain `E`, stress `S` values, and a numerical constant `C` typically calculated using sum(S ./ (E .+ 1e-10)) / length(S).

# Fields
- `E::Vector{T}`: A vector of strain values.
- `S::Vector{T}`: A vector of stress values.
- `C::T`: A computed constant based on `E` and `S` values.
"""
struct Dataset{T <: Real}
	E::Vector{T}
	S::Vector{T}
	C::T
end

"""
	struct SolveResults

Stores the results of a solving process, containing multiple properties associated with stress, strain, displacements, etc.
For all relevant fields, it stores a history of the field, i.e. one vector per iteration.

# Fields
- `N_datapoints::Int64`: Number of data points in the problem.
- `Φ::Vector{Vector{Float64}}`: Node initial coordinates or positions.
- Other fields (e.g., `e`, `E`, `s`, `S`, `u`, etc.) store different physical quantities involved in the solving process, represented as vectors.
"""
Base.@kwdef struct SolveResults
	N_datapoints::Int64
	Φ::Union{Vector{Vector{Float64}}, Vector{Float64}}
	e::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	E::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	s::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	S::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	u::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	λ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	μ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	cost::Vector{Float64} = Vector{Vector{Float64}}()
	equilibrium::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	compatibility::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
end

"""
	get_final(result::SolveResults) -> NamedTuple

Returns the final values of each field in `SolveResults`, excluding `N_datapoints` and `Φ`.

# Arguments
- `result::SolveResults`: The `SolveResults` instance to extract data from.

# Returns
- A `NamedTuple` containing the final values of all the fields except `N_datapoints` and `Φ`.
"""
function get_final(result::SolveResults)
	names = [f for f in fieldnames(SolveResults) if f ∉ ((:N_datapoints, :Φ))]
	values = [isempty(getfield(result, field)) ? [] : getfield(result, field)[end] for field in names]
	return NamedTuple{Tuple(names)}(values)
end

"""
	calc_c(E, S) -> T

Calculates the compatibility constant `C` from strain `E` and stress `S`.

# Arguments
- `E::Vector{T}`: Vector of strain values.
- `S::Vector{T}`: Vector of stress values.

# Returns
- The computed constant `C`.
"""
function calc_c(E, S)
	filter = abs.(E) .> 1e-7
	return sum(S[filter] ./ E[filter]) / sum(filter)
end

"""
	Dataset(E, S) -> Dataset

Creates a `Dataset` instance by computing the compatibility constant `C`.

# Arguments
- `E::Vector{T}`: Vector of strain values.
- `S::Vector{T}`: Vector of stress values.

# Returns
- A `Dataset` instance.
"""
function Dataset(E, S)
	return Dataset(E, S, calc_c(E, S))
end

"""
	Base.getindex(data::Dataset, i)

Accesses the `i`-th element of the `Dataset`, returning a tuple of strain and stress values.

# Arguments
- `data::Dataset`: The dataset to index.
- `i::Int`: Index to access.

# Returns
- A tuple `(E[i], S[i])`.
"""
function Base.getindex(data::Dataset, i)
	return (data.E[i], data.S[i])
end

"""
	Base.length(data::Dataset) -> Int

Returns the number of data points in the `Dataset`.

# Arguments
- `data::Dataset`: The dataset whose length is being determined.

# Returns
- The length of the dataset.
"""
function Base.length(data::Dataset)
	return length(data.E)
end

"""
	get_points(data::Dataset) -> Iterator

Returns an iterator over pairs of strain and stress values in the `Dataset`.

# Arguments
- `data::Dataset`: The dataset to iterate through.

# Returns
- An iterator of `(E, S)` pairs.
"""
function get_points(data::Dataset)
	return zip(data.E, data.S)
end

"""
	unpack_pairs(pairs) -> Tuple

Takes in a vector of `(strain, stress)` pairs and returns unpacked vectors of strain and stress values.

# Arguments
- `pairs::Vector{Tuple}`: A vector of `(strain, stress)` pairs.

# Returns
- A tuple `(strain, stress)` with individual values.
"""
function unpack_pairs(pairs)
	strain = [p[1] for p in pairs]
	stress = [p[2] for p in pairs]
	return strain, stress
end


# TODO compare to assignLocalState and choose the best implementation
"""
	choose_closest_to(target_e, target_s, w, data::Dataset) -> Tuple

Selects the closest matching `(E, S)` pairs from `data` for given target values `target_e` and `target_s`.

# Arguments
- `target_e::Vector`: Target strain values.
- `target_s::Vector`: Target stress values.
- `w::Vector`: Weights for calculating cost.
- `data::Dataset`: The dataset containing strain and stress values.

# Returns
- A tuple containing the selected `(e_values, s_values)` and the total cost.
"""
function choose_closest_to(target_e, target_s, w, data::Dataset)
	best_pairs = []  # To store the best (e, s) pairs for each element
	total_cost = 0.0
	# Iterate through each dataset in data
	for (i, (e, s)) in enumerate(zip(target_e, target_s))
		costs = collect((w[i] * (1 / 2 * data.C * (E - e)^2 + 1 / (2 * data.C) * (S - s)^2) for (E, S) in get_points(data)))
		cost, idx = findmin(costs)
		push!(best_pairs, data[idx])
		total_cost += cost
	end

	e_values, s_values = unpack_pairs(best_pairs)
	return (e_values, s_values), total_cost
end


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



function NewtonRaphsonStep(
	previous_sol::AbstractArray,
	data_star::AbstractArray,
	node_vector::AbstractArray,
	num_ele::Int,
	ndofs::AbstractArray,
	costFunc_constant::Float64,
	bar_distF::Vector{Float64},
	cross_section_area::Float64,
	alpha::Float64,
	constrained_dofs::AbstractArray;
	numQuadPts::Int = 2,
)

	# assembly
	rhs = assembleEquilibriumResidual(
		previous_sol,
		data_star,
		node_vector,
		num_ele,
		ndofs,
		costFunc_constant,
		bar_distF,
		cross_section_area,
		alpha,
		numQuadPts,
	)

	J = assembleLinearizedSystemMatrix(previous_sol, node_vector, num_ele, ndofs, costFunc_constant, cross_section_area, alpha; numQuadPts = numQuadPts)

	# enforcing boundary conditions    
	ids = collect(1:size(J, 1))
	deleteat!(ids, constrained_dofs)
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



