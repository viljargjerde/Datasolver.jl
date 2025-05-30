# import JuMP: Model, @variable, @constraint, set_attribute, @expression, @objective, backend, solve_time, set_optimizer_attribute, set_start_value, optimize!, set_silent, termination_status, MOI, value, objective_value, QuadExpr
using Gurobi
using JuMP

using MathOptInterface
const MOI = MathOptInterface

function NLP_solver(problem, dataset; use_L1_norm = true, use_data_bounds = true, random_init_data = false, verbose = false, timelimit = nothing, worklimit = nothing, write_model_file = false, parameter_file = nothing, log_intermediate_results = false)
	numDataPts = length(dataset)
	results = SolveResults(N_datapoints = numDataPts, Φ = problem.node_vector)
	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "MIPGap", 0.0)
	# set_optimizer_attribute(model, "NumericFocus", 3)
	n = problem.num_node
	m = problem.num_ele
	quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)
	N_mats, dN_mats = constructBasisFunctionMatrixLinearLagrange(1, quad_pts)
	fixed_dofs = problem.constrained_dofs[begin:length(problem.constrained_dofs)÷2]
	free_dofs = setdiff(1:n, fixed_dofs)


	if use_data_bounds
		min_s = minimum(dataset.S)
		max_s = maximum(dataset.S)
		min_e = minimum(dataset.E)
		max_e = maximum(dataset.E)
		@variable(model, min_s <= sbar[1:m] <= max_s)
		@variable(model, min_e <= ebar[1:m] <= max_e)
	else
		@variable(model, sbar[1:m])
		@variable(model, ebar[1:m])
	end
	@variable(model, uhat[1:n])


	equilibrium_eq = zeros(QuadExpr, n)
	compatibility_eq = zeros(QuadExpr, m)

	for i in 1:m
		xi0, xi1 = problem.node_vector[i:i+1]
		J4int = norm(xi1 - xi0) / 2
		J4deriv = norm(xi1 - xi0) / 2
		active_dofs_u = [i, i + 1]

		for (N_matrix, dN_mat, quad_pt, quad_weight) in zip(N_mats, dN_mats, quad_pts, quad_weights)
			integration_factor = problem.area * quad_weight * J4int
			dN_matrix = dN_mat / J4deriv
			duh = (dN_matrix*uhat[active_dofs_u])[1]

			equilibrium_eq[active_dofs_u] += (integration_factor * sbar[i] * dN_matrix' * (1.0 + problem.alpha * duh) - quad_weight * J4int * N_matrix' * problem.force((1 - quad_pt) / 2 * norm(xi0) + (1 + quad_pt) / 2 * norm(xi1)))
			compatibility_eq[i] += integration_factor * (duh + problem.alpha / 2 * duh^2 - ebar[i])
		end
	end


	@constraint(model, equilibrium_eq[free_dofs] .== 0)
	@constraint(model, compatibility_eq .== 0)

	#########
	@variable(model, choosen_ES[1:m, 1:length(dataset)], Bin)



	if !random_init_data
		s = get_initialization_s(problem)
		for ele_i in eachindex(s)
			# choose E, S pair where S is closest to s
			best_idx = argmin((abs(dataset.S[j] - s[ele_i]) for j in eachindex(dataset.S)))
			for j in eachindex(dataset.S)
				set_start_value(choosen_ES[ele_i, j], j == best_idx)
			end
		end
	end


	# * For each element exactly one datapoint has to be choosen:
	for i in 1:m
		@constraint(model, sum(choosen_ES[i, :]) == 1)
	end

	# * Fix the dofs
	for i in fixed_dofs
		@constraint(model, uhat[i] == 0)
	end

	# Define expressions for E and S to simplify later constraints
	@expression(model, E[i = 1:m], sum(dataset.E[j] * choosen_ES[i, j] for j in 1:length(dataset)))
	@expression(model, S[i = 1:m], sum(dataset.S[j] * choosen_ES[i, j] for j in 1:length(dataset)))


	if use_L1_norm
		# Define auxiliary variables for linearized absolute differences
		@variable(model, z_e[1:m] >= 0)  # Represents |e[i] - E[i]|
		@variable(model, z_s[1:m] >= 0)  # Represents |s[i] - S[i]|

		# Constraints to enforce |e[i] - E[i]| <= z_e[i] and similarly for s and S
		for i in 1:m
			@constraint(model, ebar[i] - E[i] <= z_e[i])
			@constraint(model, E[i] - ebar[i] <= z_e[i])
			@constraint(model, sbar[i] - S[i] <= z_s[i])
			@constraint(model, S[i] - sbar[i] <= z_s[i])
		end

		# Objective function minimizes the weighted sum of z_e and z_s
		@objective(model, Min, _integrateCostfunction_NLP_L1(z_e, z_s, dataset.C, problem))
	else
		@objective(model, Min, integrateCostfunction(ebar, sbar, E, S, dataset.C, problem))
		# @objective(model, Min, sum(dataset.C * (ebar[i] - E[i])^2 + 1 / dataset.C * (sbar[i] - S[i])^2 for i in 1:m))
	end
	if !verbose
		set_silent(model)
	end
	if !isnothing(timelimit)
		set_time_limit_sec(model, timelimit)
	end

	if !isnothing(worklimit)
		set_optimizer_attribute(model, "WorkLimit", worklimit)
	end

	if write_model_file
		write_to_file(model, "model.mps")
	end

	if !isnothing(parameter_file)
		open(parameter_file, "r") do f
			# Read the parameter file and set parameters in the model
			for line in eachline(f)
				if !startswith(line, "#")
					# Split the line into key and value
					key, value = split(line)
					set_optimizer_attribute(model, strip(key), parse(Int, strip(value)))
				end
			end
		end
	end

	if log_intermediate_results
		solvetime_log = Float64[]
		objective_log = Float64[]
		function results_callback(cb_data, cb_where::Cint)
			if cb_where == Gurobi.GRB_CB_MIP
				objbst = Ref{Cdouble}()
				time = Ref{Cdouble}()

				# Retrieve current best objective value
				Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_MIP_OBJBST, objbst)
				# Retrieve runtime
				Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_RUNTIME, time)
				if length(objective_log) == 0 || objective_log[end] > objbst[]
					push!(objective_log, objbst[])
					push!(solvetime_log, time[])
				end
			end
			return
		end

		# Set the callback function in the model
		MOI.set(model, Gurobi.CallbackFunction(), results_callback)
	end
	optimize!(model)


	@show solve_time(model)
	# Access Gurobi internal model
	backend = JuMP.backend(model)
	grb = backend.optimizer

	# Get work units
	work_units = MOI.get(grb, Gurobi.ModelAttribute("Work"))

	if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.TIME_LIMIT || termination_status(model) == MOI.OTHER_LIMIT
		# Get the values 
		e_values = value.(ebar)
		s_values = value.(sbar)
		u_values = value.(uhat)
		equilibrium_values = value.(equilibrium_eq[free_dofs])
		compatibility_values = value.(compatibility_eq)

		E_values = value.(E)
		S_values = value.(S)
	else
		println("The model did not solve to optimality.")
	end
	push!(results.u, u_values)
	push!(results.e, e_values)
	push!(results.s, s_values)
	push!(results.E, E_values)
	push!(results.S, S_values)
	push!(results.equilibrium, equilibrium_values)
	push!(results.compatibility, compatibility_values)
	push!(results.solvetime, work_units)
	push!(results.solvetime, solve_time(model))
	if termination_status(model) == MOI.TIME_LIMIT || termination_status(model) == MOI.OTHER_LIMIT
		push!(results.solvetime, 0.0)
		printstyled("The model reached the time limit.\n"; color = :red)
	else
		push!(results.solvetime, 1.0)
	end
	push!(results.cost, integrateCostfunction(e_values, s_values, E_values, S_values, dataset.C, problem, L2 = !use_L1_norm))

	if log_intermediate_results
		append!(results.solvetime, solvetime_log[2:end])
		append!(results.cost, objective_log[2:end])

	end
	return results
end


function _integrateCostfunction_NLP_L1(abs_ediff::AbstractArray, abs_sdiff::AbstractArray, costFunc_constant::Float64, problem::Barproblem)

	# quad points in default interval [-1,1]
	_, quad_weights = GaussLegendreQuadRule(numQuadPts = problem.num_quad_pts)
	# integration
	costFunc_global = 0.0
	for i in 1:problem.num_ele      # loop over element
		# jacobian for the integration
		xi0, xi1 = problem.node_vector[i:i+1]
		J4int = norm(xi1 - xi0) / 2
		costFunc_global += 0.5 * (sqrt(costFunc_constant) * abs_ediff[i] + 1 / sqrt(costFunc_constant) * abs_sdiff[i]) * sum(quad_weights) * J4int * problem.area
	end

	return costFunc_global
end
