import LinearAlgebra.I
import GaussQuadrature.legendre
using SparseArrays

"""
Args:
	constraints: A vector of tuples, where each tuple contains a node ID and the degree of freedom (DOF) to be constrained. Both indices start from 1.
		Example: (3,2) means the 3rd node is constrained in the y-direction.
	num_elements: The total number of elements in the system.
	dimensions: The number of spatial dimensions (1, 2, or 3).

Returns:
	A vector of constrained degrees of freedom.
"""
function get_constrained_dofs(constraints::Vector{Tuple{Int64, Int64}}, num_ele::Int64, dims::Int64)

	ndof_u = (num_ele + 1) * dims
	ndof_e = ndof_s = ndof_mu = num_ele

	lambda_offset = ndof_u + ndof_e + ndof_s + ndof_mu

	fixed_dofs = zeros(Int64, length(constraints) * 2) # Both for u and also λ

	for i in eachindex(constraints)
		(node, dof) = constraints[i]
		@assert dof <= dims "Invalid degree of freedom fixed, check dimensions"
		local_dof = (node - 1) * dims + dof
		fixed_dofs[i] = local_dof
		fixed_dofs[i+length(constraints)] = lambda_offset + local_dof

	end
	return fixed_dofs

end


function constructBasisFunctionMatrixLinearLagrange(
	dims::Int,
	quad_pts::AbstractVector;
	interval::AbstractArray = [-1, 1],
)
	N0, N1 = linearLagrangePolynomials(interval)
	dN0, dN1 = compute1stDeriv4linearLagrangePolynomials(interval)

	I_mat = I(dims)  # Identity matrix

	N_mats = Vector{Matrix{Float64}}(undef, length(quad_pts))
	dN_mats = Vector{Matrix{Float64}}(undef, length(quad_pts))

	for (i, quad_pt) in enumerate(quad_pts)
		N_mats[i] = hcat(N0(quad_pt) * I_mat, N1(quad_pt) * I_mat)
		dN_mats[i] = hcat(dN0(quad_pt) * I_mat, dN1(quad_pt) * I_mat)
	end

	return N_mats, dN_mats
end

# function constructBasisFunctionMatrixLinearLagrange(dims::Int64; interval::AbstractArray = [-1, 1])
# 	N0, N1 = linearLagrangePolynomials(interval)
# 	dN0, dN1 = compute1stDeriv4linearLagrangePolynomials(interval)

# 	I_mat = I(dims)  # Identity matrix


# 	return x -> hcat([N0(x) * I_mat, N1(x) * I_mat]...), x -> (hcat([dN0(x) * I_mat, dN1(x) * I_mat]...))
# end


# linear Lagrange polynomial on an interval [x0,x1]
function linearLagrangePolynomials(interval::AbstractArray)
	x0, x1 = interval
	L0 = x -> (x - x1) / (x0 - x1)
	L1 = x -> (x - x0) / (x1 - x0)

	return L0, L1
end


function compute1stDeriv4linearLagrangePolynomials(interval::AbstractArray)
	x0, x1 = interval
	dL0 = x -> 1 / (x0 - x1)
	dL1 = x -> 1 / (x1 - x0)

	return dL0, dL1
end



# Gauss-Legendre quadrature rule

function GaussLegendreQuadRule(; interval::AbstractArray = [-1, 1], numQuadPts::Int = 2)
	# Gauss-Legendre quadrature rule on [-1,1]
	x, w = legendre(numQuadPts)

	# map to an arbitrary interval
	x0, x1 = interval
	m, n = (x1 - x0) / 2, (x1 + x0) / 2

	return x .* m .+ n, w .* m
end





function assembleEquilibriumResidual(
	x,
	E,
	S,
	costFunc_constant,
	problem,
	load_alpha,
)
	dims = problem.dims
	# quad points in default interval [-1,1]
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)

	# basis function matrix evaluated in master element [-1,1]
	N_mats, dN_mats = constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)

	# extract variable fields from solution vector of the previous iteration
	ndofs = get_ndofs(problem)
	ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs
	rhs = zeros(sum(ndofs))
	r_u = 1:ndof_u
	r_e = last(r_u)+1:last(r_u)+ndof_e
	r_s = last(r_e)+1:last(r_e)+ndof_s
	r_mu = last(r_s)+1:last(r_s)+ndof_mu
	r_lambda = last(r_mu)+1:last(r_mu)+ndof_lambda

	uhat      = @view x[r_u]
	ebar      = @view x[r_e]
	sbar      = @view x[r_s]
	mubar     = @view x[r_mu]
	lambdahat = @view x[r_lambda]

	# prepare the difference between material data
	e_diff = @views ebar - E
	s_diff = @views sbar - S

	# alloccation blocks of rhs
	rhs_b1 = @view rhs[r_u]
	rhs_b2 = @view rhs[r_e]
	rhs_b3 = @view rhs[r_s]
	rhs_b4 = @view rhs[r_mu]
	rhs_b5 = @view rhs[r_lambda]

	alpha = problem.alpha

	# assembly routine
	@views for cc_ele ∈ 1:problem.num_ele      # loop over elements    
		active_dofs_u = active_dofs_lambda = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
		active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele

		# jacobian for the integration
		xi0, xi1 = problem.node_vector[cc_ele:cc_ele+1]
		J4int = norm(xi1 - xi0) / 2

		# jacobian for derivative
		J4deriv = norm(xi1 - xi0) / 2

		eh = ebar[active_dofs_e]
		sh = sbar[active_dofs_s]
		muh = mubar[active_dofs_mu]

		for (N_matrix, dN_mat, quad_pt, quad_weight) in zip(N_mats, dN_mats, quad_pts, quad_weights)

			integration_factor = problem.area * quad_weight * J4int

			dN_matrix = dN_mat / J4deriv

			dPhih = dN_matrix * [xi0; xi1]
			duh = dN_matrix * uhat[active_dofs_u]

			dlambdah = dN_matrix * lambdahat[active_dofs_lambda]
			e_uh = duh ⋅ dPhih + alpha / 2 * duh ⋅ duh
			PBh = (dPhih + alpha * duh)

			# integrated blocks of the rhs
			rhs_b1[active_dofs_u] += integration_factor * (-alpha * dN_matrix' * sh * dlambdah -
														   dN_matrix' * PBh * muh)

			rhs_b2[active_dofs_e] += integration_factor * (muh - costFunc_constant * e_diff[active_dofs_e])
			rhs_b3[active_dofs_s] += integration_factor * (-PBh' * dlambdah - s_diff[active_dofs_s] / costFunc_constant)

			rhs_b4[active_dofs_mu] += -integration_factor * (e_uh - eh)
			rhs_b5[active_dofs_lambda] += N_matrix' * quad_weight * J4int * problem.force(quad_pt) * load_alpha -
										  dN_matrix' * integration_factor * PBh * sh
		end
	end

	return rhs
end



function assembleLinearizedSystemMatrix(x, problem::Barproblem, costFunc_constant::Float64)
	dims = problem.dims
	# quad points in default interval [-1,1]
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)

	# basis function matrix evaluated in master element [-1,1]
	_, dN_mats = constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)

	ndofs = get_ndofs(problem)
	ndof_total = sum(ndofs)
	ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs

	r_u = 1:ndof_u
	r_e = last(r_u)+1:last(r_u)+ndof_e
	r_s = last(r_e)+1:last(r_e)+ndof_s
	r_mu = last(r_s)+1:last(r_s)+ndof_mu
	r_lambda = last(r_mu)+1:last(r_mu)+ndof_lambda

	uhat      = @view x[r_u]
	sbar      = @view x[r_s]
	mubar     = @view x[r_mu]
	lambdahat = @view x[r_lambda]

	# Allocate the full sparse matrix
	J = spzeros(ndof_total, ndof_total)

	# Compute index ranges for each variable block
	r_u      = 1:ndof_u
	r_e      = last(r_u)+1:last(r_u)+ndof_e
	r_s      = last(r_e)+1:last(r_e)+ndof_s
	r_mu     = last(r_s)+1:last(r_s)+ndof_mu
	r_lambda = last(r_mu)+1:last(r_mu)+ndof_lambda

	# Define views for the upper triangle blocks (explicitly computed)
	J11 = @view J[r_u, r_u]
	J13 = @view J[r_u, r_s]
	J14 = @view J[r_u, r_mu]
	J15 = @view J[r_u, r_lambda]

	J22 = @view J[r_e, r_e]
	J24 = @view J[r_e, r_mu]

	J33 = @view J[r_s, r_s]
	J35 = @view J[r_s, r_lambda]



	# Define views for the symmetric lower triangle blocks (transposes of the above)
	J31 = @view J[r_s, r_u]
	J41 = @view J[r_mu, r_u]
	J42 = @view J[r_mu, r_e]
	J51 = @view J[r_lambda, r_u]
	J53 = @view J[r_lambda, r_s]


	alpha = problem.alpha

	# assembly routine
	@views for cc_ele ∈ 1:problem.num_ele      # loop over elements
		active_dofs_u = active_dofs_lambda = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
		active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele
		xi0, xi1 = problem.node_vector[cc_ele:cc_ele+1]
		# jacobian for the integration
		J4int = norm(xi1 - xi0) / 2
		# jacobian for derivative
		J4deriv = norm(xi1 - xi0) / 2

		sh = sbar[active_dofs_s]
		muh = mubar[active_dofs_mu]

		for (dN_mat, quad_weight) in zip(dN_mats, quad_weights)

			dN_matrix = dN_mat / J4deriv
			integration_factor = problem.area * quad_weight * J4int

			duh = dN_matrix * uhat[active_dofs_u]

			dPhih = dN_matrix * [xi0; xi1]
			PBh = (dPhih + alpha * duh)

			dlambdah = dN_matrix * lambdahat[active_dofs_lambda]

			# integrated blocks of the system matrix
			J11[active_dofs_u, active_dofs_u] += alpha * integration_factor * muh * dN_matrix' * dN_matrix

			J13[active_dofs_u, active_dofs_s] += alpha * integration_factor * dN_matrix' * dlambdah

			J14[active_dofs_u, active_dofs_mu] += integration_factor * dN_matrix' * (dPhih + alpha * duh)

			J15[active_dofs_u, active_dofs_lambda] += alpha * integration_factor * sh * dN_matrix' * dN_matrix

			J22[active_dofs_e, active_dofs_e] += integration_factor * costFunc_constant

			J24[active_dofs_e, active_dofs_mu] += -integration_factor

			J33[active_dofs_s, active_dofs_s] += integration_factor / costFunc_constant

			J35[active_dofs_s, active_dofs_lambda] += (integration_factor*PBh'*dN_matrix)[:]
		end
	end


	J31 .= J13'
	J41 .= J14'
	J42 .= J24'
	J51 .= J15'
	J53 .= J35'


	return J
end



