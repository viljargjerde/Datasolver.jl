export NewtonRaphsonStep
import LinearAlgebra.I


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

	fixed_dofs = Vector()

	for (node, dof) in constraints
		@assert dof <= dims "Invalid degree of freedom fixed, check dimensions"
		local_dof = (node - 1) * dims + dof
		push!(fixed_dofs, local_dof)
		push!(fixed_dofs, lambda_offset + local_dof)

	end
	sort!(fixed_dofs)
	return fixed_dofs

end

function reconstruct_vector(vec, constrained_dofs)
	original_length = length(vec) + length(constrained_dofs)
	reconstructed = zeros(original_length)  # Initialize with zeros
	j = 1  # Index for vec
	k = 1  # Index for constrained_dofs
	num_constraints = length(constrained_dofs)

	for i in 1:original_length
		if k <= num_constraints && i == constrained_dofs[k]
			# Insert a zero at constrained indices
			k += 1
		else
			# Copy elements from vec
			reconstructed[i] = vec[j]
			j += 1
		end
	end

	return reconstructed
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




export linearLagrangePolynomials, compute1stDeriv4linearLagrangePolynomials, constantFunctions

export constructBasisFunctionMatrixLinearLagrange, constructBasisFunctionMatrixConstantFuncs

export GaussLegendreQuadRule

import GaussQuadrature.legendre


function constructBasisFunctionMatrixLinearLagrange(dims::Int64; interval::AbstractArray = [-1, 1])
	N0, N1 = linearLagrangePolynomials(interval)
	dN0, dN1 = compute1stDeriv4linearLagrangePolynomials(interval)

	I_mat = I(dims)  # Identity matrix


	return x -> sparse(hcat([N0(x) * I_mat, N1(x) * I_mat]...)), x -> (sparse(hcat([dN0(x) * I_mat, dN1(x) * I_mat]...)))
end

# TODO Check if this can be removed
function constructBasisFunctionMatrixConstantFuncs(; evalPts::AbstractArray = [-1, 1])
	R = constantFunctions()

	R_matrix = Matrix((R.(evalPts))')

	return sparse(R_matrix)
end


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


function constantFunctions()
	L = x -> 1.0

	return L
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


export assembleLinearSystemMatrix, assembleRhsLinearBar
export assembleEquilibriumResidual, assembleLinearizedSystemMatrix



function assembleEquilibriumResidual(
	x::AbstractArray,
	data_star::AbstractArray,
	node_vector::AbstractArray,
	num_ele::Int,
	ndofs::AbstractArray,
	costFunc_constant::Float64,
	bar_distF::Vector{Float64},
	cross_section_area::Float64,
	alpha::Float64,
	numQuadPts::Int,
)
	dims = length(node_vector[1])
	# quad points in default interval [-1,1]
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)

	# basis function matrix evaluated in master element [-1,1]
	N_func, dN_func = constructBasisFunctionMatrixLinearLagrange(dims;)
	# R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts=quad_pts)

	# extract variable fields from solution vector of the previous iteration
	ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs

	uhat = x[1:ndof_u]
	ebar = x[ndof_u+1:ndof_u+ndof_e]
	sbar = x[ndof_u+ndof_e+1:ndof_u+ndof_e+ndof_s]
	mubar = x[ndof_u+ndof_e+ndof_s+1:ndof_u+ndof_e+ndof_s+ndof_mu]
	lambdahat = x[ndof_u+ndof_e+ndof_s+ndof_mu+1:end]

	# prepare the difference between material data
	e_diff = ebar - data_star[:, 1]
	s_diff = sbar - data_star[:, 2]

	# alloccation blocks of rhs
	rhs_b1 = spzeros(ndof_u)
	rhs_b5 = spzeros(ndof_lambda)

	rhs_b2 = spzeros(ndof_e)
	rhs_b3 = spzeros(ndof_s)
	rhs_b4 = spzeros(ndof_mu)

	# assembly routine
	for cc_ele ∈ 1:num_ele      # loop over elements    
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)
			active_dofs_u = active_dofs_lambda = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
			active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele

			# jacobian for the integration
			xi0, xi1 = node_vector[cc_ele:cc_ele+1]
			J4int = norm(xi1 - xi0) / 2


			# jacobian for derivative
			J4deriv = norm(xi1 - xi0) / 2

			# interpolate discrete solution from the previous iteration

			integration_factor = cross_section_area * quad_weight * J4int
			N_matrix = N_func(quad_pt)
			# @show J4deriv

			dN_matrix = dN_func(quad_pt) / J4deriv
			dPhih = dN_matrix * [xi0; xi1]
			duh = dN_matrix * uhat[active_dofs_u]

			eh = ebar[active_dofs_e]
			sh = sbar[active_dofs_s]
			muh = mubar[active_dofs_mu]

			dlambdah = dN_matrix * lambdahat[active_dofs_lambda]
			e_uh = duh ⋅ dPhih + alpha / 2 * duh ⋅ duh
			PBh = (dPhih + alpha * duh)

			# integrated blocks of the rhs
			rhs_b1[active_dofs_u] += integration_factor * (-alpha * dN_matrix' * sh * dlambdah -
														   dN_matrix' * PBh * muh)

			rhs_b2[active_dofs_e] += integration_factor * (muh - costFunc_constant * e_diff[active_dofs_e])
			rhs_b3[active_dofs_s] += integration_factor * (-PBh' * dlambdah - s_diff[active_dofs_s] / costFunc_constant)

			rhs_b4[active_dofs_mu] += -integration_factor * (e_uh - eh)
			rhs_b5[active_dofs_lambda] += N_matrix' * quad_weight * J4int * bar_distF -
										  dN_matrix' * integration_factor * PBh * sh
		end
	end

	# global rhs
	rhs = [rhs_b1; rhs_b2; rhs_b3; rhs_b4; rhs_b5]

	return rhs
end



function assembleLinearizedSystemMatrix(x::AbstractArray, node_vector::AbstractArray, num_ele::Int, ndofs::AbstractArray, costFunc_constant::Float64, cross_section_area::Float64, alpha::Float64; numQuadPts::Int = 2)
	dims = size(node_vector[1])[1]
	# quad points in default interval [-1,1]
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = numQuadPts)

	# basis function matrix evaluated in master element [-1,1]
	_, dN_func = constructBasisFunctionMatrixLinearLagrange(dims)


	# extract variable fields from solution vector of the previous iteration
	ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs

	uhat = x[1:ndof_u]
	sbar = x[ndof_u+ndof_e+1:ndof_u+ndof_e+ndof_s]
	mubar = x[ndof_u+ndof_e+ndof_s+1:ndof_u+ndof_e+ndof_s+ndof_mu]
	# mubar = rand(ndof_mu)
	lambdahat = x[ndof_u+ndof_e+ndof_s+ndof_mu+1:end]
	# lambdahat = rand(ndof_lambda)

	# alloccation blocks of the system matrix
	J11, J12, J13, J14, J15 = spzeros(ndof_u, ndof_u), spzeros(ndof_u, ndof_e), spzeros(ndof_u, ndof_s), spzeros(ndof_u, ndof_mu), spzeros(ndof_u, ndof_lambda)

	J22, J23, J24, J25 = spzeros(ndof_e, ndof_e), spzeros(ndof_e, ndof_s), spzeros(ndof_e, ndof_mu), spzeros(ndof_e, ndof_lambda)

	J33, J34, J35 = spzeros(ndof_s, ndof_s), spzeros(ndof_s, ndof_mu), spzeros(ndof_s, ndof_lambda)

	J44, J45 = spzeros(ndof_mu, ndof_mu), spzeros(ndof_mu, ndof_lambda)

	J55 = spzeros(ndof_lambda, ndof_lambda)

	# assembly routine
	for cc_ele ∈ 1:num_ele      # loop over elements
		for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)
			active_dofs_u = active_dofs_lambda = collect((cc_ele-1)*dims+1:(cc_ele+1)*dims)
			active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele

			# jacobian for the integration
			xi0, xi1 = node_vector[cc_ele:cc_ele+1]
			J4int = norm(xi1 - xi0) / 2

			# jacobian for derivative
			J4deriv = norm(xi1 - xi0) / 2
			# @show J4deriv

			dN_matrix = dN_func(quad_pt) / J4deriv
			integration_factor = cross_section_area * quad_weight * J4int

			# interpolate discrete solution from the previous iteration
			duh = dN_matrix * uhat[active_dofs_u]

			sh = sbar[active_dofs_s]
			muh = mubar[active_dofs_mu]
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


	# global system matrix
	J = [
		J11 J12 J13 J14 J15;
		J12' J22 J23 J24 J25;
		J13' J23' J33 J34 J35;
		J14' J24' J34' J44 J45;
		J15' J25' J35' J45' J55
	]

	return J
end

