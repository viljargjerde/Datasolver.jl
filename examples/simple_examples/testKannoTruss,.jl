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
