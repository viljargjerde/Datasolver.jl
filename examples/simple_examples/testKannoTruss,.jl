using Datasolver, Revise, LinearAlgebra, Test

###### TEST TRUSS STRUCTURE ######

@testset "Recovery of 1d fixed-free linear straight bar" begin
    bar_L = 2*3.6
    A = 2000/1e6        # [m²]
    bar_E = 1.622e+07   # [Pa] 

    alpha = 0.0         # 0: linear strain     1: nonlinear strain

    λ = 1               # load factor
    F = 400.0*λ         # [N]

    uxRef = F*bar_L/A/bar_E
    eRef = F/A/bar_E
    sRef = F/A

    node_vector = [
        [0,         0],
        [bar_L/2,   0],
        [bar_L,     0]
    ]

    constrained_dofs = [
        (1, 1),
        (1, 2)
    ]

    connections = [
        (1, 2),
        (2, 3)
    ]

    force = zeros(2 * length(node_vector))
    force[end-1] = F


    problem = TrussProblem(
        A,
        force,
        connections,
        alpha,
        constrained_dofs,
        node_vector = node_vector,
        num_quad_pts = 2,
    )

    num_ele = length(node_vector)-1
    num_node = length(node_vector)
    dims = length(node_vector[1])

    ndof_u = ndof_lambda = num_node * dims
    ndof_e = ndof_s = ndof_mu = num_ele
    ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
    ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]

    free_dofs = collect(1:ndof_tot)
    deleteat!(free_dofs, problem.constrained_dofs)

    x = zeros(ndof_tot)
    E = [eRef;eRef]
    S = [sRef;sRef]

    # assembly
    rhs = Datasolver.assembleEquilibriumResidual(
        x,
        E,
        S,
        bar_E,
        problem,
    );

    J = Datasolver.assembleLinearizedSystemMatrix(x, problem, bar_E);

    # enforcing boundary conditions    
    J_free = J[free_dofs, free_dofs]
    rhs_free = rhs[free_dofs]

    # solving
    Delta_x = zero(x)
    Delta_x[free_dofs] = Datasolver.qr(J_free) \ rhs_free

    # check residual
    r::Float64 = norm(rhs_free - J_free * Delta_x[free_dofs])
    @test r <= 1e-12

    # check linear solution
    xRef = zeros(ndof_tot)
    xRef[1:2:ndof_u] = (F/A/bar_E) .* [node_vector[i][1] for i in 1:problem.num_node]
    xRef[ndof_u+1:ndof_u+ndof_e] .= eRef
    xRef[ndof_u+ndof_e+1:ndof_u+ndof_e+ndof_s] .= sRef

    @test isapprox(Delta_x, xRef)
end

@testset "2D linear truss" begin
    A = 2000/1e6        # [m²]
    bar_E = 1.622e+09   # [Pa] 
    
    alpha = 0.0         # 0: linear strain     1: nonlinear strain
    NR_max_iter = 10
    NR_tol = 2e-10
    
    λ = 1               # load factor
    F = -400.0*λ         # [N]
    
    # reference solution computed with https://valdivia.staff.jade-hs.de/fachwerk_en.html
    uRef = [0;0;0;-1.072e-3;0;0;0;-6.278e-4]
    sRef = [0;0;-sqrt(2)*1e+5;2.0e+5;-sqrt(2)*1e+5]
    eRef = sRef ./ bar_E
    
    node_vector = [
        [0,     0],
        [3.6,   0],
        [7.2,   0],
        [3.6,   3.6]
    ]
    
    constrained_dofs = [
        (1, 1),
        (1, 2),
        (3, 1),
        (3, 2)
    ]
    
    connections = [
        (1, 2),
        (2, 3),
        (1, 4),
        (2, 4),
        (3, 4)
    ]
    
    force = zeros(2 * length(node_vector))
    force[4] = F
    
    
    problem = TrussProblem(
        A,
        force,
        connections,
        alpha,
        constrained_dofs,
        node_vector = node_vector,
        num_quad_pts = 2,
    )
    
    num_ele = problem.num_ele
    num_node = problem.num_node
    dims = problem.dims
    
    ndof_u = ndof_lambda = num_node * dims
    ndof_e = ndof_s = ndof_mu = num_ele
    ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
    ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]
    
    free_dofs = collect(1:ndof_tot)
    deleteat!(free_dofs, problem.constrained_dofs)
    
    x = zeros(ndof_tot)
    E = eRef
    S = sRef
    ii = 0
    
    for iter in 1:NR_max_iter
        ii += 1
        Delta_x = Datasolver.NewtonRaphsonStep(
                        x,
                        E,
                        S,
                        bar_E,
                        problem,
                        free_dofs,
                        true,
                    )
    
        # update solution
        x += Delta_x
    
        # check convergence
        if norm(Delta_x) <= NR_tol
            break
        end
    end
    
    
    xRef = zeros(ndof_tot)
    xRef[1:ndof_u+ndof_e+ndof_s] = [uRef;eRef;sRef]
    
    err = x-xRef
    
    @test ii == 2
    @test norm(err) <= 5e-7    
end
