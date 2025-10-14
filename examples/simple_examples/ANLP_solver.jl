
using LinearAlgebra, SparseArrays, Revise, Datasolver



function checkThermomechanicalConsistency(;results::SolveResults)

    # check nonnegative work of the chosen data points
    w = results.E[end] .* results.S[end]
    id1_thermoinconsistent = findall(x -> x < 0, w)

    if isempty(id1_thermoinconsistent) == false
        println("The work of chosen data points in elements $id1_thermoinconsistent is negative. These elements are thermomechanical inconsistent.")
    end

    # check nonnegative work of the final phase state
    w = results.e[end] .* results.s[end]
    id2_thermoinconsistent = findall(x -> x < 0, w)

    if isempty(id2_thermoinconsistent) == false
        println("The work of computed strain and stress in elements $id2_thermoinconsistent is negative. These elements are thermomechanical inconsistent.")
    end

    # check nonnegative product of chosen stress data and computed stress field
    w = results.s[end] .* results.S[end]
    id3_thermoinconsistent = findall(x -> x < 0, w)

    if isempty(id3_thermoinconsistent) == false
        println("The product of chosen stress data and computed stress in elements $id3_thermoinconsistent is negative. These elements are thermomechanical inconsistent.")
    end

    if isempty(id1_thermoinconsistent) && isempty(id2_thermoinconsistent) && isempty(id3_thermoinconsistent)
        println("All elements of the discretized structure are thermomechanical consistent.")
    end

end


#--------------FUNCTIONALITIES FOR THE STANDARD 3-FIELDS MIXED FORMULATION--------------------


function solveANLP(;
    initProblem::Dataproblem,
    constrained_dofs_global,
    externalForce,
    num_load_steps::Int64=1,
    loadFac::Vector{Float64}=[1.0],
    YoungModulus::Float64=1.0,
    scaleFacYoungModulus::Float64=1.0,
    NR_tol::Float64=1e-10,
    NR_max_iter::Int=100,
    qrFactorized::Bool=true)

    # allocation + pre-processing
    node_vector = initProblem.node_vector
    results = SolveResults(N_datapoints=1, Φ=node_vector)
    
    ndofs = Datasolver.get_ndofs(initProblem)
    ndof_u, ndof_e, ndof_s = ndofs[1:3]
    ndof_tot = ndof_u + ndof_e + ndof_s

    # initial guess of the solution for the 1st load step
    x = zeros(ndof_tot)
    NRiter = Int64[]

    for i = 1:num_load_steps
        λl = loadFac[i+1] * scaleFacYoungModulus

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

        # newton-raphson scheme
        x, iter, flag = NewtonRaphsonSchemeANLP(problem = problem, 
                                              currentSol = x,
                                              NR_tol = NR_tol,
                                              NR_max_iter = NR_max_iter, 
                                              YoungModulus = YoungModulus * scaleFacYoungModulus,
                                              qrFactorized = qrFactorized)
        
        # collect computed solution fields
        indices = cumsum(ndofs)

        uhat = x[1:indices[1]]
        ebar = x[indices[1]+1:indices[2]]
        sbar = x[indices[2]+1:indices[3]] ./ scaleFacYoungModulus
        
        push!(results.u, collect(uhat))
        push!(results.e, collect(ebar))
        push!(results.s, collect(sbar))
        push!(NRiter, iter)

        if flag == 1
            println("NR scheme did not converge at load step $i. Computation stopped.")
            push!(results.NRiter, NRiter)
            break
        end
    end

    push!(results.NRiter, NRiter)

    nn, ii = findmax(NRiter)
    println("Computation finished. Max NRiter = $nn at load step $ii.")

    return results
end



function NewtonRaphsonSchemeANLP(;problem::Dataproblem, 
                                currentSol::AbstractArray,
                                YoungModulus::Float64=1.0,
                                NR_tol::Float64=1e-10,
                                NR_max_iter::Int64=50,
                                qrFactorized::Bool=true)

    iter = 0
    flag = 0

    while iter <= NR_max_iter
        # solve the linearized system
        Δx = solveLinSys(problem=problem, currentSol=currentSol, YoungModulus=YoungModulus, qrFactorized=qrFactorized)

        # update the solution
        currentSol += Δx

        # check convergence
        if norm(Δx) <= NR_tol
            break
        end

        iter += 1
    end

    if iter >= NR_max_iter
        flag = 1
    end

    return currentSol, iter, flag
end



function solveLinSys(;problem::Dataproblem, 
                     currentSol::AbstractArray,
                     YoungModulus::Float64=1.0,
                     qrFactorized::Bool=true)

    # allocate unconstrained solution vector
    ndofs = Datasolver.get_ndofs(problem)
    ndof_u, ndof_e, ndof_s = ndofs[1:3]
    ndof_tot = ndof_u + ndof_e + ndof_s

    numConDofs = Int64(ceil(length(problem.constrained_dofs) / 2))
    free_dofs = collect(1:ndof_tot)
    deleteat!(free_dofs, problem.constrained_dofs[1:numConDofs])

    # assembly constrained matrix equations
    A = assembleLinSysMatrixANLP(problem, currentSol, YoungModulus)
    b = assembleResidualsANLP(problem, currentSol, YoungModulus)

    Ac = A[free_dofs,free_dofs]
    rhs = b[free_dofs]

    # solve
    Δx = zeros(ndof_tot)
    if qrFactorized
        Δx[free_dofs] = qr(Ac) \ rhs
    else
        Δx[free_dofs] = Ac \ rhs
    end

    # check residual and condition number of the system matrix
    r::Float64 = norm(rhs - Ac * Δx[free_dofs])
    κ::Float64 = cond(Matrix(Ac))

    if r > 1e-10
        println("Warning: Solution with residual $r > 1e-10")
    end
    if κ > 1e20
        # @show κ
        println("Warning: Condition number of the system matrix $κ > 1e20")
    end

    return Δx
end


function assembleResidualsANLP(problem, currentSol, YoungModulus)
    # quad points in default interval [-1,1]
    quad_pts, quad_weights = Datasolver.GaussLegendreQuadRule(numQuadPts=problem.num_quad_pts)

    # basis function matrix evaluated in master element [-1,1]
    dims = problem.dims
    N_mats, dN_mats = Datasolver.constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)

    # extract variable fields from the current solution vector
    ndofs = Datasolver.get_ndofs(problem)
    ndof_u, ndof_e, ndof_s = ndofs[1:3]
    ndof_tot = ndof_u + ndof_e + ndof_s

    r_u = 1:ndof_u
    r_e = last(r_u)+1:last(r_u)+ndof_e
    r_s = last(r_e)+1:last(r_e)+ndof_s
    
    uhat = @view currentSol[r_u]
    ebar = @view currentSol[r_e]
    sbar = @view currentSol[r_s]
    
    # alloccation blocks of rhs
    rhs = zeros(ndof_tot)    
    rhs_b1 = @view rhs[r_u]
    rhs_b2 = @view rhs[r_e]
    rhs_b3 = @view rhs[r_s]

    # assembly routine
    α = problem.alpha

    @views for cc_ele ∈ 1:problem.num_ele      # loop over elements
        # indices of active dofs
        ele_a, ele_b = problem.connections[cc_ele]
        active_dofs_u = vcat((ele_a-1) * dims+1 : ele_a * dims,
                             (ele_b-1) * dims+1 : ele_b * dims)
        active_dofs_e = active_dofs_s = cc_ele

        # jacobian for the integration
        xi0 = problem.node_vector[ele_a]
        xi1 = problem.node_vector[ele_b]
        J4int = norm(xi1 - xi0) / 2

        # jacobian for derivative
        J4deriv = norm(xi1 - xi0) / 2

        eh = ebar[active_dofs_e]
        sh = sbar[active_dofs_s]
        
        for (N_matrix, dN_mat, quad_pt, quad_weight) in zip(N_mats, dN_mats, quad_pts, quad_weights)            # loop over quadrature points

            integration_factor = problem.area * quad_weight * J4int            

            dN_matrix = dN_mat / J4deriv

            dPhih = dN_matrix * [xi0; xi1]
            duh = dN_matrix * uhat[active_dofs_u]

            e_uh = dot(duh, dPhih) + α/2 .* dot(duh, duh)
            PBh = (dPhih + α .* duh)

            # integrated blocks of the rhs
            rhs_b1[active_dofs_u] += - (dN_matrix' * PBh) .* (integration_factor * sh)

            if problem.force isa Function
                x_quad = (1 - quad_pt) / 2 * norm(xi0) + (1 + quad_pt) / 2 * norm(xi1)
                rhs_b1[active_dofs_u] += (quad_weight * J4int) .* (N_matrix' * problem.force(x_quad))
            end

            rhs_b2[active_dofs_e] += integration_factor * (sh - YoungModulus * eh)
            rhs_b3[active_dofs_s] += integration_factor * (eh - e_uh)
        end             # end loop over quadrature points
    end                 # end loop over elements

    if (problem.force isa Function) == false
        rhs[r_u] += problem.force
    end

    return rhs
end


function assembleLinSysMatrixANLP(problem, currentSol, YoungModulus)
    # quad points in default interval [-1,1]
    quad_pts, quad_weights = Datasolver.GaussLegendreQuadRule(numQuadPts=problem.num_quad_pts)

    # basis function matrix evaluated in master element [-1,1]
    dims = problem.dims
    N_mats, dN_mats = Datasolver.constructBasisFunctionMatrixLinearLagrange(dims, quad_pts)

    # extract variable fields from the current solution vector
    ndofs = Datasolver.get_ndofs(problem)
    ndof_u, ndof_e, ndof_s = ndofs[1:3]
    ndof_tot = ndof_u + ndof_e + ndof_s

    r_u = 1:ndof_u
    r_e = last(r_u)+1:last(r_u)+ndof_e
    r_s = last(r_e)+1:last(r_e)+ndof_s
    
    uhat = @view currentSol[r_u]
    sbar = @view currentSol[r_s]

    # allocation
    K = spzeros(ndof_tot,ndof_tot)

    K11 = @view K[r_u, r_u]
    K13 = @view K[r_u, r_s]

    K22 = @view K[r_e, r_e]
    K23 = @view K[r_e, r_s]

    K31 = @view K[r_s, r_u]
    K32 = @view K[r_s, r_e]

    # assembly routine
    α = problem.alpha

    @views for cc_ele ∈ 1:problem.num_ele      # loop over elements
        # indices of active dofs
        ele_a, ele_b = problem.connections[cc_ele]
        active_dofs_u = vcat((ele_a-1) * dims+1 : ele_a * dims,
                             (ele_b-1) * dims+1 : ele_b * dims)
        active_dofs_e = active_dofs_s = cc_ele

        # jacobian for the integration
        xi0 = problem.node_vector[ele_a]
        xi1 = problem.node_vector[ele_b]
        J4int = norm(xi1 - xi0) / 2

        # jacobian for derivative
        J4deriv = norm(xi1 - xi0) / 2

        sh = sbar[active_dofs_s]

        for (dN_mat, quad_weight) in zip(dN_mats, quad_weights)            # loop over quadrature points
            integration_factor = problem.area * quad_weight * J4int            

            dN_matrix = dN_mat / J4deriv

            dPhih = dN_matrix * [xi0; xi1]
            duh = dN_matrix * uhat[active_dofs_u]
            PBh = (dPhih + α .* duh)

            # integrated blocks
            K11[active_dofs_u, active_dofs_u] += (integration_factor * α * sh) .* (dN_matrix' * dN_matrix)

            K13[active_dofs_u, active_dofs_s] += integration_factor .* (dN_matrix' * PBh)

            K22[active_dofs_e, active_dofs_e] += integration_factor * YoungModulus

            K23[active_dofs_e, active_dofs_s] += -integration_factor

            K31[active_dofs_s, active_dofs_u] += integration_factor .* ((PBh' * dN_matrix)[:])

            K32[active_dofs_s, active_dofs_e] += -integration_factor
        end                                    # end loop over quadrature points
    end                                        # end loop over elements

    return K
end

