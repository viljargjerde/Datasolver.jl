import JuMP: Model, @variable, @constraint, @expression, @objective, optimize!, termination_status, MOI, value, objective_value
using Gurobi


# function LP_solver(connections, Φ, A, data, f, fixed_dofs; verbose=false)
#     results = SolveResults(N_datapoints=length(data), Φ=Φ)
#     B = create_B_matrix(connections, Φ)
#     B = remove_dofs(B, fixed_dofs)
#     w = calc_w(connections, Φ, A)
#     m, n = size(B) # Elements, nodes
#     N_data = length(data)
#     if length(f) != n
#         f = remove_dofs(f, fixed_dofs)
#     end
#     model = Model(Gurobi.Optimizer)

#     @variable(model, u[1:n])
#     @variable(model, s[1:m])
#     @variable(model, e[1:m])
#     @variable(model, choosen_ES[1:m, 1:N_data], Bin)

#     # * For each element exactly one datapoint has to be choosen:
#     for i in 1:m
#         @constraint(model, sum(choosen_ES[i, :]) == 1)
#         @constraint(model, sum(choosen_ES[i, :]) == 1)
#     end

#     # * Balance and compatibility
#     @constraint(model, B' * (s .* w) .== f)
#     @constraint(model, e .== B * u)


#     @expression(model, E[i=1:m], sum(data.E[j] * choosen_ES[i, j] for j in 1:N_data))
#     @expression(model, S[i=1:m], sum(data.S[j] * choosen_ES[i, j] for j in 1:N_data))



#     # Define auxiliary variables for linearized absolute differences
#     @variable(model, z_e[1:m] >= 0)  # Represents |e[i] - E[i]|
#     @variable(model, z_s[1:m] >= 0)  # Represents |s[i] - S[i]|

#     # Constraints to enforce |e[i] - E[i]| <= z_e[i] and similarly for s and S
#     for i in 1:m
#         @constraint(model, e[i] - E[i] <= z_e[i])
#         @constraint(model, E[i] - e[i] <= z_e[i])
#         @constraint(model, s[i] - S[i] <= z_s[i])
#         @constraint(model, S[i] - s[i] <= z_s[i])
#     end

#     # Now set up a new objective function that minimizes the weighted sum of z_e and z_s
#     @objective(model, Min, sum(sqrt(w[i]) * (sqrt(data.C) * z_e[i] + sqrt((1 / data.C)) * z_s[i]) for i in 1:m))



#     optimize!(model)


#     if termination_status(model) == MOI.OPTIMAL
#         # Get the values of e and s
#         e_values = value.(e)
#         s_values = value.(s)
#         u_values = value.(u)

#         # Get the selected E and S values for each index
#         E_values = [value(E[i]) for i in 1:m]
#         S_values = [value(S[i]) for i in 1:m]

#         balance = B' * (s_values .* w) - f
#         compatibility = e_values - B * u_values

#         # Display the results
#     else
#         println("The model did not solve to optimality.")
#     end
#     push!(results.e, e_values)
#     push!(results.s, s_values)
#     push!(results.E, E_values)
#     push!(results.S, S_values)
#     push!(results.u, u_values)
#     push!(results.balance, balance)
#     push!(results.compatibility, compatibility)
#     push!(results.cost, objective_value(model))
#     return results
# end

function nonlin_LP_solver(node_vector, data, cross_section_area, bar_distF, fixed_dofs; alpha=1.0, use_L1_norm=true)
    numDataPts = length(data)
    results = SolveResults(N_datapoints=numDataPts, Φ=node_vector)
    model = Model(Gurobi.Optimizer)
    n = length(node_vector)
    m = n - 1
    quad_pts, quad_weights = DataDrivenNonlinearBar.GaussLegendreQuadRule(numQuadPts=2)
    N_func, dN_func = DataDrivenNonlinearBar.constructBasisFunctionMatrixLinearLagrange(1)
    J4ints = [norm(node_vector[i+1] - node_vector[i]) / 2 for i in 1:m]
    J4derivs = [norm(node_vector[i+1] - node_vector[i]) / 2 for i in 1:m]
    @variable(model, uhat[1:n])
    @variable(model, sbar[1:m])
    @variable(model, ebar[1:m])
    @variable(model, mubar[1:m])
    @variable(model, lambdahat[1:n])

    internal_force = (sum((cross_section_area * quad_weight * J4ints[i]) * sbar[i] * (1.0 + alpha * (dN_func(quad_pt)/J4derivs[i]*uhat[i:i+1])[1]) for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)) for i in 1:m)

    external_force = (sum((quad_weight * J4ints[i]) * (N_func(quad_pt)*bar_distF)[1] for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)) for i in 1:m)

    compatibility_eq = (sum(
        (cross_section_area * quad_weight * J4ints[i]) * ((dN_func(quad_pt)/J4derivs[i]*uhat[i:i+1])[1] + alpha / 2 * (dN_func(quad_pt)/J4derivs[i]*uhat[i:i+1])[1]^2 - ebar[i]) for (quad_pt, quad_weight) in zip(quad_pts, quad_weights)) for i in 1:m
    )


    @constraint(model, internal_force .== external_force)
    @constraint(model, compatibility_eq .== 0)

    #########
    @variable(model, choosen_ES[1:m, 1:length(data)], Bin)

    # * For each element exactly one datapoint has to be choosen:
    for i in 1:m
        @constraint(model, sum(choosen_ES[i, :]) == 1)
        @constraint(model, sum(choosen_ES[i, :]) == 1)
    end

    # * Fix the dofs
    for i in fixed_dofs[begin:length(fixed_dofs)÷2]
        @constraint(model, uhat[i] == 0)
        @constraint(model, lambdahat[i] == 0)
    end

    # Define expressions for E and S to simplify later constraints
    @expression(model, E[i=1:m], sum(data.E[j] * choosen_ES[i, j] for j in 1:length(data)))
    @expression(model, S[i=1:m], sum(data.S[j] * choosen_ES[i, j] for j in 1:length(data)))


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
        @objective(model, Min, sum((sqrt(data.C) * z_e[i] + sqrt(1 / data.C) * z_s[i]) for i in 1:m))
    else
        @objective(model, Min, sum(data.C * (ebar[i] - E[i])^2 + 1 / data.C * (sbar[i] - S[i])^2) for i in 1:m)
    end


    optimize!(model)


    if termination_status(model) == MOI.OPTIMAL
        # Get the values of e and s
        e_values = value.(ebar)
        s_values = value.(sbar)
        u_values = value.(uhat)
        λ_values = value.(lambdahat)
        μ_values = value.(mubar)

        # Get the selected E and S values for each index
        E_values = [value(E[i]) for i in 1:m]
        S_values = [value(S[i]) for i in 1:m]


        # Display the results
    else
        println("The model did not solve to optimality.")
    end
    push!(results.e, e_values)
    push!(results.s, s_values)
    push!(results.E, E_values)
    push!(results.S, S_values)
    push!(results.u, u_values)
    push!(results.λ, λ_values)
    push!(results.μ, μ_values)
    push!(results.cost, objective_value(model))
    return results
end
