


using LinearAlgebra


costFunc_ele(e, s, C) = 0.5 * (C * e^2 + 1 / C * s^2);
costFunc_ele_L1(e, s, C) = 0.5 * (sqrt(C) * abs(e) + sqrt(1 / C) * abs(s));


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


function find_closest_idx(S::Vector{Float64}, s::Vector{Float64})
    idxs = zeros(Int64, length(s))
    for i in eachindex(s)
        idxs[i] = argmin((abs(S[j] - s[i]) for j in eachindex(S)))
    end
    return idxs
end

function directSolverNonLinearBar(
    problem::Dataproblem,
    dataset::Dataset;
    init_indices=nothing,
    random_init_data::Bool=false,
    DD_max_iter::Int=100,
    NR_tol::Float64=1e-10,
    NR_max_iter::Int=100,
    verbose::Bool=false,
)
    start_time = time()

    ## initialize e_star and s_star
    numDataPts = length(dataset)
    node_vector = problem.node_vector
    num_ele = problem.num_ele
    results = SolveResults(N_datapoints=numDataPts, Φ=node_vector)
    num_node = length(node_vector)
    dims = length(node_vector[1])
    ndof_u = ndof_lambda = num_node * dims
    ndof_e = ndof_s = ndof_mu = num_ele

    ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
    ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]
    free_dofs = collect(1:ndof_tot)
    deleteat!(free_dofs, problem.constrained_dofs)

    data_idxs_old = Int64[]
    if init_indices !== nothing
        E = dataset.E[init_indices]
        S = dataset.S[init_indices]
        data_idxs_old = deepcopy(init_indices)
    elseif random_init_data
        init_data_id = rand(1:numDataPts, num_ele)
        E = dataset.E[init_data_id]
        S = dataset.S[init_data_id]
        data_idxs_old = deepcopy(init_data_id)
    else
        s = get_initialization_s(problem)
        best_idxs = find_closest_idx(dataset.S, s)
        S = dataset.S[best_idxs]
        E = dataset.E[best_idxs]
        data_idxs_old = deepcopy(best_idxs)
    end
    # iterative data-driven direct solver
    x = zeros(ndof_tot)
    dd_iter = 0
    NRiter = Int64[]
    start_solvetime = time()
    while dd_iter <= DD_max_iter

        # newton-raphson scheme
        cc_iter = 0
        for iter in 1:NR_max_iter
            cc_iter += 1
            Delta_x = NewtonRaphsonStep(
                x,
                E,
                S,
                dataset.C,
                problem,
                free_dofs,
                verbose,
            )

            # update solution
            x += Delta_x

            # check convergence
            if norm(Delta_x) <= NR_tol
                break
            end

            if iter == NR_max_iter && verbose
                println("NR did not converge")
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
        data_idxs = assignLocalState(dataset, ebar, sbar)

        new_E = dataset.E[data_idxs]
        new_S = dataset.S[data_idxs]
        curr_cost = integrateCostfunction(ebar, sbar, E, S, dataset.C, problem)

        converged = data_idxs_old == data_idxs
        dd_iter += 1

        # overwrite local state
        E = new_E
        S = new_S
        equilibrium = equilibrium_eq(uhat, sbar, problem)
        compat = compatibility_eq(uhat, ebar, problem)
        
        push!(results.u, collect(uhat))
        push!(results.e, collect(ebar))
        push!(results.s, collect(sbar))
        push!(results.λ, [norm(λ[i:i+dims-1]) for i in 1:dims:length(λ)])
        push!(results.μ, collect(μ))
        push!(results.E, collect(E))
        push!(results.S, collect(S))

        push!(results.data_idx, data_idxs)
        push!(results.cost, curr_cost)
        push!(results.equilibrium, equilibrium)
        push!(results.compatibility, compat)

        push!(NRiter, cc_iter)

        if converged
            end_time = time()
            push!(results.solvetime, end_time - start_time)
            push!(results.solvetime, end_time - start_solvetime)
            @assert norm(equilibrium) < NR_tol "norm(equilibrium) = $(norm(equilibrium))"
            @assert norm(compat) < NR_tol "norm(compatibility) = $(norm(compat))"
            break
        else
            data_idxs_old = deepcopy(data_idxs)
        end
    end

    push!(results.NRiter, NRiter)
    push!(results.ADMiter, dd_iter)

    nn = maximum(results.NRiter)
    println("Computation takes $dd_iter ADM iters and up to $nn NR iters.")

    @assert argmin(results.cost) == length(results.cost) "The last result is not the best one"
    return results
end


function greedyLocalSearchSolverNonLinearBar(
    problem::Dataproblem,
    dataset::Dataset;
    random_init_data::Bool=false,
    DD_max_iter::Int=100,
    NR_tol::Float64=1e-10,
    NR_max_iter::Int=100,
    verbose::Bool=false,
    search_iters::Int=100,
    cache_ADM::Bool=true,
)
    start_time = time()
    result = SolveResults(N_datapoints=length(dataset), Φ=problem.node_vector)
    ADM_cache = cache_ADM ? Set{Vector{Int64}}() : nothing

    first_result = directSolverNonLinearBar(problem, dataset;
        random_init_data=random_init_data,
        DD_max_iter=DD_max_iter,
        NR_tol=NR_tol,
        NR_max_iter=NR_max_iter,
        verbose=verbose,
    )
    push_final_result!(result, first_result)
    push!(result.solvetime, time() - start_time)
    if cache_ADM
        for d_idx in first_result.data_idx
            push!(ADM_cache, d_idx)
        end
    end

    search_iter = 1
    while search_iter <= search_iters
        diffs = costFunc_ele.(result.E[end] - result.e[end], result.S[end] - result.s[end], dataset.C)
        sorted_idx = sortperm(diffs, rev=true)  # biggest first


        for j in sorted_idx
            search_iter += 1
            trial_data_idxs = copy(result.data_idx[end])

            # Try finding the closest index for this specific element
            local_diffs = costFunc_ele.(dataset.E .- result.e[end][j], dataset.S .- result.s[end][j], dataset.C)
            min_idx1, min_idx2 = find_two_smallest_indices(local_diffs)
            if trial_data_idxs[j] == min_idx1
                trial_data_idxs[j] = min_idx2
            else
                trial_data_idxs[j] = min_idx1
            end
            if cache_ADM && in(trial_data_idxs, ADM_cache)
                if verbose
                    println("Skip this trial, already computed")
                end
                continue  # skip if already computed
            end
            trial_result = directSolverNonLinearBar(
                problem,
                dataset;
                init_indices=trial_data_idxs,
                DD_max_iter=DD_max_iter,
                NR_tol=NR_tol,
                NR_max_iter=NR_max_iter,
                verbose=verbose,
            )
            if cache_ADM
                for d_idx in trial_result.data_idx
                    push!(ADM_cache, d_idx)
                end
            end
            if trial_result.cost[end] < result.cost[end]
                # accept move
                push_final_result!(result, trial_result)
                push!(result.solvetime, time() - start_time)
                break  # restart from the top
            end

            if j == sorted_idx[end]
                # no improving move found
                search_iter = search_iters + 1
            end
            if search_iter > search_iters
                break
            end
        end
    end
    end_time = time()
    push!(result.solvetime, end_time - start_time)
    return result
end


function find_two_smallest_indices(vec::Vector{<:Real})
    if vec[1] < vec[2]
        min1, min2 = vec[1], vec[2]
        idx1, idx2 = 1, 2
    else
        min1, min2 = vec[2], vec[1]
        idx1, idx2 = 2, 1
    end

    for i in 3:length(vec)
        if vec[i] < min1
            min2, idx2 = min1, idx1
            min1, idx1 = vec[i], i
        elseif vec[i] < min2
            min2, idx2 = vec[i], i
        end
    end

    return (idx1, idx2)
end


function equilibrium_eq(uhat, sbar, problem::Dataproblem)
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=problem.num_quad_pts)
    dims = problem.dims
    N_mats, dN_mats = constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)
    equilibrium = zeros((problem.num_node) * dims)
    alpha = problem.alpha
    for cc_ele ∈ 1:problem.num_ele      # loop over elements  
		ele_a, ele_b = problem.connections[cc_ele]  
		active_dofs_u = vcat((ele_a-1)*dims+1 : ele_a * dims,
                     (ele_b-1) * dims+1 : ele_b * dims)
        for (N_matrix, dN_mat, quad_pt, quad_weight) in zip(N_mats, dN_mats, quad_pts, quad_weights)
			ele_a, ele_b = problem.connections[cc_ele]
            active_dofs_s = cc_ele
            sh = sbar[active_dofs_s]

            # jacobian for the integration
            xi0 = problem.node_vector[ele_a]
            xi1 = problem.node_vector[ele_b]
            J4int = norm(xi1 - xi0) / 2

            # jacobian for derivative
            J4deriv = norm(xi1 - xi0) / 2
            dN_matrix = dN_mat / J4deriv
            duh = dN_matrix * uhat[active_dofs_u]
            dPhih = dN_matrix * [xi0; xi1]

            PBh = (dPhih + alpha * duh)
            integration_factor = problem.area * quad_weight * J4int

            if problem.force isa Function
				equilibrium[active_dofs_u] += N_matrix' * (quad_weight * J4int * problem.force((1 - quad_pt) / 2 * norm(xi0) + (1 + quad_pt) / 2 * norm(xi1))) -
                                          (dN_matrix') * (integration_factor) * (PBh * sh)
			else
				equilibrium[active_dofs_u] += - dN_matrix' * integration_factor * (PBh * sh)
			end
        end
    end

    if (problem.force isa Function) == false
        equilibrium += problem.force
    end

    idxs = collect(1:length(equilibrium))
    deleteat!(idxs, problem.constrained_dofs[begin:length(problem.constrained_dofs)÷2])
    equilibrium[idxs] # Remove constrained dofs
end


function compatibility_eq(uhat, ebar, problem::Dataproblem)
    dims = problem.dims
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=problem.num_quad_pts)
    _, dN_mats = constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)
    compatibility = zeros(problem.num_ele)
    alpha = problem.alpha

    for cc_ele ∈ 1:problem.num_ele      # loop over elements    
        for (dN_mat, quad_weight) in zip(dN_mats, quad_weights)
			ele_a, ele_b = problem.connections[cc_ele]
			xi0 = problem.node_vector[ele_a]
        	xi1 = problem.node_vector[ele_b]
			active_dofs_u = vcat((ele_a-1)*dims+1 : ele_a * dims,
                     (ele_b-1) * dims+1 : ele_b * dims)
            # jacobian for the integration
            J4int = norm(xi1 - xi0) / 2
            # jacobian for derivative
            J4deriv = norm(xi1 - xi0) / 2
            dN_matrix = dN_mat / J4deriv

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

    # # allocation
    indices = zeros(Int64, length(ebar))
    # find the closest data point to the local state
    Threads.@threads for i ∈ eachindex(indices)

        distances = (costFunc_ele(dataset.E[j] - ebar[i], dataset.S[j] - sbar[i], dataset.C) for j in 1:length(dataset))
        indices[i] = argmin(distances)
    end

    return indices
end



function get_initialization_s(problem::Dataproblem)
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=problem.num_quad_pts)
    dims = problem.dims
    N_mats, dN_mats = constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)
    ndof_u, _, ndof_s, _, _ = get_ndofs(problem)

    # Construct global matrices A sbar = b
    A = zeros(ndof_u, ndof_s)
    b = zeros(ndof_u)
    for cc_ele ∈ 1:problem.num_ele      # loop over elements 
		ele_a, ele_b = problem.connections[cc_ele]   
        active_dofs_u = vcat((ele_a-1)*dims+1 : ele_a * dims,
                     (ele_b-1) * dims+1 : ele_b * dims)
      
		active_dofs_s = cc_ele
        # jacobian for the integration
        xi0 = problem.node_vector[ele_a]
		xi1= problem.node_vector[ele_b] 
        J4int = norm(xi1 - xi0) / 2

        # jacobian for derivative
        J4deriv = norm(xi1 - xi0) / 2
        for (N_matrix, dN_mat, quad_pt, quad_weight) in zip(N_mats, dN_mats, quad_pts, quad_weights)
            dN_matrix = dN_mat / J4deriv
            dPhih = dN_matrix * [xi0; xi1]

            integration_factor = problem.area * quad_weight * J4int
            A[active_dofs_u, active_dofs_s] += integration_factor * dN_matrix' * dPhih
            if problem.force isa Function
                b[active_dofs_u] += N_matrix' * quad_weight * J4int * problem.force((1 - quad_pt) / 2 * norm(xi0) + (1 + quad_pt) / 2 * norm(xi1))
            end
        end

    end

    if (problem.force isa Function) == false
        b += problem.force
    end

    idxs = collect(1:length(b))
    deleteat!(idxs, problem.constrained_dofs[begin:length(problem.constrained_dofs)÷2]) # constrained_dofs include lambda, but here we only care about u, which is the first half
    @views A[idxs, :] \ b[idxs]
end


function integrateCostfunction(e::AbstractArray, s::AbstractArray, E::AbstractArray, S::AbstractArray, costFunc_constant::Float64, problem::Dataproblem; L2=true)

    # quad points in default interval [-1,1]
    _, quad_weights = GaussLegendreQuadRule(numQuadPts=problem.num_quad_pts)

    # integration
    costFunc_global = 0.0

    for i in 1:problem.num_ele      # loop over element
        # jacobian for the integration
		ele_a, ele_b = problem.connections[i]
		xi0 = problem.node_vector[ele_a]
        xi1 = problem.node_vector[ele_b]
        J4int = norm(xi1 - xi0) / 2
        if L2
            costFunc_global += costFunc_ele(e[i] - E[i], s[i] - S[i], costFunc_constant) * sum(quad_weights) * J4int * problem.area

        else
            costFunc_global += costFunc_ele_L1(e[i] - E[i], s[i] - S[i], costFunc_constant) * sum(quad_weights) * J4int * problem.area
        end
    end

    return costFunc_global
end



function NewtonRaphsonStep(
    x::AbstractArray,
    E::AbstractArray,
    S::AbstractArray,
    costFunc_constant::Float64,
    problem::Dataproblem,
    free_dofs::AbstractArray,
    verbose::Bool,
    QRfactorized::Bool=true
)

    # assembly

    rhs = assembleEquilibriumResidual(
        x,
        E,
        S,
        costFunc_constant,
        problem,
    )

    J = assembleLinearizedSystemMatrix(x, problem, costFunc_constant)

    # enforcing boundary conditions    
    J_free = J[free_dofs, free_dofs]
    rhs_free = rhs[free_dofs]
    # solving
    Delta_x = zero(x)
    if QRfactorized
        Delta_x[free_dofs] = qr(J_free) \ rhs_free
    else
        Delta_x[free_dofs] = J_free \ rhs_free
    end


    if verbose
        # check residual and condition number of J
        r::Float64 = norm(rhs_free - J_free * Delta_x[free_dofs])
        κ::Float64 = cond(Matrix(J_free))

        if r > 1e-10
            println("Warning: Solution with residual $r > 1e-10")
        end
        if κ > 1e20
            # @show κ
            println("Warning: Condition number of the system matrix $κ > 1e20")
        end
    end
    return Delta_x
end




function directSolverNonLinearBarA(;
        initProblem::Dataproblem,
        constrained_dofs_global,
        externalForce,
        dataset::Dataset,
        scaleFactorDataConst::Float64=1.0,
        num_load_steps::Int64=1,
        loadFac::Vector{Float64}=[1.0],
        init_indices=nothing,
        random_init_data::Bool=false,
        solution_guess1stload=nothing,
        DD_max_iter::Int=100,
        NR_tol::Float64=1e-10,
        NR_max_iter::Int=100,
        verbose::Bool=false,
        QRfactorized::Bool=true
    )

    # allocation
    numDataPts = length(dataset)
    node_vector = initProblem.node_vector
    results = SolveResults(N_datapoints=numDataPts, Φ=node_vector)

    num_ele = initProblem.num_ele
    ndofs = Datasolver.get_ndofs(initProblem)
    ndof_tot = sum(ndofs)

    free_dofs = collect(1:ndof_tot)
    deleteat!(free_dofs, initProblem.constrained_dofs)

    # initial guess of the solution for the 1st load step
    if solution_guess1stload !== nothing
        x = solution_guess1stload
    else
        x = zeros(ndof_tot)
    end
    
    E = Float64[]
    S = Float64[]
    data_idxs_old = Int64[]

    for i = 1:num_load_steps
        λl = loadFac[i+1] * scaleFactorDataConst

        global E
        global S
        global data_idxs_old

        if externalForce isa Function
            Ffunc = x -> externalForce(x,λl)
            Fnodal = [0]
        else
            Ffunc = x -> 0.0
            Fnodal = externalForce .* λl
        end

        problem = TrussProblem(
            initProblem.area,
            Fnodal,
            initProblem.connections,
            initProblem.alpha,
            constrained_dofs_global,
            node_vector = node_vector,
            num_quad_pts = initProblem.num_quad_pts,
            force_func = Ffunc
        )

        # initialize e_star and s_star
        if i == 1            
            if init_indices !== nothing
                E = dataset.E[init_indices]
                S = dataset.S[init_indices]
                data_idxs_old = deepcopy(init_indices)
            elseif random_init_data
                init_data_id = rand(1:numDataPts, num_ele)
                E = dataset.E[init_data_id]
                S = dataset.S[init_data_id]
                data_idxs_old = deepcopy(init_data_id)
            else
                s = Datasolver.get_initialization_s(problem)
                best_idxs = Datasolver.find_closest_idx(dataset.S, s)
                S = dataset.S[best_idxs]
                E = dataset.E[best_idxs]
                data_idxs_old = deepcopy(best_idxs)
            end
        else
            # using e_star and s_star of the previous load step
            E = results.E[end]
            S = results.S[end]
            data_idxs_old = results.data_idx[end]
        end

        # iterative data-driven direct solver    
        x, results = directSolverNonLinearBarB!(
            problem = problem,
            results = results,
            currentSol = x,
            activeDofsIds = free_dofs,
            dataE = E,
            dataS = S,
            dataset = dataset,
            data_idxs_current = data_idxs_old,
            scaleFactorDataConst = scaleFactorDataConst,
            DD_max_iter=DD_max_iter,
            NR_max_iter=NR_max_iter,
            NR_tol=NR_tol,
            verbose=verbose,
            QRfactorized=QRfactorized
            )

        dd_iter = results.ADMiter[i]
        nn = maximum(results.NRiter)
        println("Computation takes $dd_iter ADM iters and up to $nn NR iters at load step $i.")
    end

    return results
end


function directSolverNonLinearBarB!(;
    problem::Dataproblem,
    results::SolveResults,
    currentSol::AbstractArray,
    activeDofsIds::AbstractArray,
    dataE::AbstractArray,
    dataS::AbstractArray,
    dataset::Dataset,
    data_idxs_current::AbstractArray,
    scaleFactorDataConst::Float64=1.0,
    DD_max_iter::Int=100,
    NR_max_iter::Int=50,
    NR_tol::Float64=1e-10,
    verbose::Bool=false,
    QRfactorized::Bool=true
    )

    x = currentSol      # currentSol is updated when x is updated!

    ndofs = Datasolver.get_ndofs(problem)
    indices = cumsum(ndofs)

    dims = problem.dims

    # iterative data-driven direct solver    
    dd_iter = 0
    NRiter = Int64[]

    while dd_iter <= DD_max_iter
        # newton-raphson scheme
        cc_iter = 0
        for iter in 1:NR_max_iter
            cc_iter += 1
            Delta_x = Datasolver.NewtonRaphsonStep(
                x,
                dataE,
                dataS,
                dataset.C,
                problem,
                activeDofsIds,
                verbose,
                QRfactorized
            )
        
            # update solution
            x += Delta_x
        
            # check convergence
            if norm(Delta_x) <= NR_tol
                break
            end
        
            if iter == NR_max_iter && verbose
                println("NR did not converge")
                break
            end
        end

        # collect computed ebar and sbar: extract variables from x using computed indices
        uhat = x[1:indices[1]]
        ebar = x[indices[1]+1:indices[2]]
        sbar = x[indices[2]+1:indices[3]]
        μ = x[indices[3]+1:indices[4]]
        λ = x[indices[4]+1:end]

        ## local state assignment
        data_idxs = Datasolver.assignLocalState(dataset, ebar, sbar)

        new_E = dataset.E[data_idxs]
        new_S = dataset.S[data_idxs]
        curr_cost = Datasolver.integrateCostfunction(ebar, sbar, dataE, dataS, dataset.C, problem)
        
        converged = data_idxs_current == data_idxs
        dd_iter += 1

        # overwrite local state
        dataE = new_E
        dataS = new_S

        equilibrium = Datasolver.equilibrium_eq(uhat, sbar, problem)
        compat = Datasolver.compatibility_eq(uhat, ebar, problem)

        push!(results.u, collect(uhat))
        push!(results.e, collect(ebar))
        push!(results.s, collect(sbar) ./ scaleFactorDataConst)
        push!(results.λ, collect(λ))
        push!(results.μ, collect(μ) ./ scaleFactorDataConst)
        push!(results.E, collect(dataE))
        push!(results.S, collect(dataS))
    
        push!(results.data_idx, data_idxs)
        push!(results.cost, curr_cost)
        push!(results.equilibrium, equilibrium)
        push!(results.compatibility, compat)
    
        push!(NRiter, cc_iter)
            
        if converged
            @assert norm(equilibrium) < NR_tol "norm(equilibrium) = $(norm(equilibrium))"
            @assert norm(compat) < NR_tol "norm(compatibility) = $(norm(compat))"
            break
        else
            data_idxs_current = deepcopy(data_idxs)
        end
    end

    push!(results.NRiter, NRiter)
    push!(results.ADMiter, dd_iter)

    return currentSol, results 
end



function greedyLocalSearchSolverNonLinearBarA(;
    initProblem::Dataproblem,
    constrained_dofs_global,
    externalForce,
    dataset::Dataset,
    scaleFactorDataConst::Float64=1.0,
    num_load_steps::Int64=1,
    loadFac::Vector{Float64}=[1.0],
    init_indices=nothing,
    random_init_data::Bool=false,
    solution_guess1stload=nothing,
    DD_max_iter::Int=100,
    NR_tol::Float64=1e-10,
    NR_max_iter::Int=100,
    verbose::Bool=false,
    search_iters::Int=100,
    cache_ADM::Bool=true,
    QRfactorized::Bool=true
    )
    
    # allocation
    numDataPts = length(dataset)
    node_vector = initProblem.node_vector
    results = SolveResults(N_datapoints=numDataPts, Φ=node_vector)

    num_ele = initProblem.num_ele
    ndofs = Datasolver.get_ndofs(initProblem)
    ndof_tot = sum(ndofs)

    free_dofs = collect(1:ndof_tot)
    deleteat!(free_dofs, initProblem.constrained_dofs)

    # initial guess of the solution for the 1st load step
    if solution_guess1stload !== nothing
        x = solution_guess1stload
    else
        x = zeros(ndof_tot)
    end
    
    E = Float64[]
    S = Float64[]
    data_idxs_old = Int64[]

    start_time = time()
    
    for i = 1:num_load_steps
        λl = loadFac[i+1] * scaleFactorDataConst

        global E
        global S
        global data_idxs_old

        if externalForce isa Function
            Ffunc = x -> externalForce(x,λl)
            Fnodal = [0]
        else
            Ffunc = x -> 0.0
            Fnodal = externalForce .* λl
        end

        problem = TrussProblem(
            initProblem.area,
            Fnodal,
            initProblem.connections,
            initProblem.alpha,
            constrained_dofs_global,
            node_vector = node_vector,
            num_quad_pts = initProblem.num_quad_pts,
            force_func = Ffunc
        )

        # initialize e_star and s_star for the 1st result of first load step
        if i == 1            
            if init_indices !== nothing
                E = dataset.E[init_indices]
                S = dataset.S[init_indices]
                data_idxs_old = deepcopy(init_indices)
            elseif random_init_data
                init_data_id = rand(1:numDataPts, num_ele)
                E = dataset.E[init_data_id]
                S = dataset.S[init_data_id]
                data_idxs_old = deepcopy(init_data_id)
            else
                s = Datasolver.get_initialization_s(problem)
                best_idxs = Datasolver.find_closest_idx(dataset.S, s)
                S = dataset.S[best_idxs]
                E = dataset.E[best_idxs]
                data_idxs_old = deepcopy(best_idxs)
            end
        else
            # using e_star and s_star of the previous load step
            E = results.E[end]
            S = results.S[end]
            data_idxs_old = results.data_idx[end]
        end

        # GO-ADM solver
        println("Load step $i:")

        x, results = greedyLocalSearchSolverNonLinearBarB(
            problem = problem,
            results = results,
            currentSol = x,
            activeDofsIds = free_dofs,
            dataE = E,
            dataS = S,
            dataset = dataset,
            data_idxs_current = data_idxs_old,
            scaleFactorDataConst = scaleFactorDataConst,
            DD_max_iter=DD_max_iter,
            NR_max_iter=NR_max_iter,
            NR_tol=NR_tol,
            verbose=verbose,
            QRfactorized=QRfactorized,
            search_iters = search_iters,
            cache_ADM = cache_ADM
        )

        dd_iter = results.ADMiter[i]
        nn = maximum(results.NRiter)
        println("   Computation takes up to $dd_iter ADM iters and $nn NR iters.")
    end

    end_time = time()
    push!(results.solvetime, end_time - start_time)

    return results
end


function greedyLocalSearchSolverNonLinearBarB(;
    problem::Dataproblem,
    results::SolveResults,
    currentSol::AbstractArray,
    activeDofsIds::AbstractArray,
    dataE::AbstractArray,
    dataS::AbstractArray,
    dataset::Dataset,
    data_idxs_current::AbstractArray,
    scaleFactorDataConst::Float64 = 1.0,
    DD_max_iter::Int=100,
    NR_max_iter::Int=50,
    NR_tol::Float64=1e-10,
    verbose::Bool=false,
    QRfactorized::Bool=true,
    search_iters::Int=100,
    cache_ADM::Bool=true
    )

    # allocation for GO-ADM results at the current load step
    result_i = SolveResults(N_datapoints=length(dataset), Φ=problem.node_vector)

    first_result = SolveResults(N_datapoints=length(dataset), Φ=problem.node_vector)

    ADM_cache = cache_ADM ? Set{Vector{Int64}}() : nothing

    start_time = time()

    # copy the current solution and phase state of the current load step
    # these will be updated after finishing solving
    x = deepcopy(currentSol)
    E, S = deepcopy(dataE), deepcopy(dataS)
    data_idxs = deepcopy(data_idxs_current)

    # first result at the current load step
    x, first_result = directSolverNonLinearBarB!(
                        problem=problem,
                        results=first_result,
                        currentSol=x,
                        activeDofsIds=activeDofsIds,
                        dataE=E,
                        dataS=S,
                        dataset=dataset,
                        data_idxs_current=data_idxs,
                        scaleFactorDataConst = scaleFactorDataConst,
                        DD_max_iter=DD_max_iter,
                        NR_max_iter=NR_max_iter,
                        NR_tol=NR_tol,
                        verbose=verbose,
                        QRfactorized=QRfactorized
                        )
    
    push_final_result!(result_i, first_result)
    push!(result_i.solvetime, time() - start_time)
    if cache_ADM
        for d_idx in first_result.data_idx
            push!(ADM_cache, d_idx)
        end
    end

    # "greedy" search loop
    search_iter = 1
    while search_iter <= search_iters
        diffs = costFunc_ele.(result_i.E[end] - result_i.e[end], result_i.S[end] - result_i.s[end] .* scaleFactorDataConst, dataset.C)
        sorted_idx = sortperm(diffs, rev=true)  # biggest first

        for j in sorted_idx     # loop over elements (starting with max cost function value)
            search_iter += 1
            trial_data_idxs = copy(result_i.data_idx[end])

            # Try finding the closest index for this specific element
            local_diffs = costFunc_ele.(dataset.E .- result_i.e[end][j], dataset.S .- result_i.s[end][j] * scaleFactorDataConst, dataset.C)
            min_idx1, min_idx2 = find_two_smallest_indices(local_diffs)

            if trial_data_idxs[j] == min_idx1
                trial_data_idxs[j] = min_idx2
            else
                trial_data_idxs[j] = min_idx1
            end

            if cache_ADM && in(trial_data_idxs, ADM_cache)
                if verbose
                    println("Skip this trial, already computed")
                end
                continue  # skip if already computed
            end

            # recompute with new data as initial data
            x = deepcopy(currentSol)
            E = dataset.E[trial_data_idxs]
            S = dataset.S[trial_data_idxs]
            data_idxs = deepcopy(trial_data_idxs)

            trial_result = SolveResults(N_datapoints=length(dataset), Φ=problem.node_vector)

            x, trial_result = directSolverNonLinearBarB!(
                        problem=problem,
                        results=trial_result,
                        currentSol=x,
                        activeDofsIds=activeDofsIds,
                        dataE=E,
                        dataS=S,
                        dataset=dataset,
                        data_idxs_current=data_idxs,
                        scaleFactorDataConst = scaleFactorDataConst,
                        DD_max_iter=DD_max_iter,
                        NR_max_iter=NR_max_iter,
                        NR_tol=NR_tol,
                        verbose=verbose,
                        QRfactorized=QRfactorized
                        )

            if cache_ADM
                for d_idx in trial_result.data_idx
                    push!(ADM_cache, d_idx)
                end
            end

            # comparing cost function
            if trial_result.cost[end] < result_i.cost[end]
                # accept move
                push_final_result!(result_i, trial_result)
                push!(result_i.solvetime, time() - start_time)
                break  # restart from the top
            end

            # if no improvement found
            if j == sorted_idx[end]
                # no improving move found
                search_iter = search_iters + 1
            end
            if search_iter > search_iters
                break
            end
        end         # end loop over elements (starting with max cost function value)
    end             # end "greedy" search loop

    end_time = time()
    push!(result_i.solvetime, end_time - start_time)

    # add the GO-ADM optimized results of the current load step output
    push_final_result!(results, result_i)
    currentSol = deepcopy(x)

    println("   $search_iter searches for GO-ADM,")    

    return currentSol, results
end
