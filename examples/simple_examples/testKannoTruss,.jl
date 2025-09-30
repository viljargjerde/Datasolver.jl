using Datasolver, Revise, LinearAlgebra, Test, Plots, LaTeXStrings

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

@testset "Recovery 1D nonlinear bar with manifactured solution" begin    
    bar_L = Float64(π)
    A = 2000/1e6        # [m²]
    bar_E = 1.622e+03   # [Pa] 
    
    α = 1.0         # 0: linear strain     1: nonlinear strain
    
    ne = 8
    
    numDataPts = 512
    
    NR_max_iter = 50
    NR_tol = 1e-8
    
    # manifactured solution and nonlinear force function
    β = 0.1*π
    
    uRef(x,β) = β * sin(π*x/bar_L)
    
    eRef(x,α,β) = β * π / bar_L * cos(π*x/bar_L) * ( 1 + 0.5*α * β * π / bar_L * cos(π*x/bar_L) )
    
    sRef(x,α,β) = bar_E * eRef(x,α,β)
    
    force_func(α, β, x, λ) = λ .*
        [bar_E * A * ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_L))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_L)) / bar_L))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_L)) / bar_L) + 2) / bar_L^(2)); 0]
    
    
    # mesh
    h = bar_L/ne
    node_vector = [ [(i-1)*h, 0] for i in 1:ne+1 ]
    
    constrained_dofs = [
        (1, 1),
        (1, 2),
        (ne+1,1),
        (ne+1,2)
    ]
    
    connections = [ (i, i+1) for i in 1:ne ]
    
    # strain limit = max of eRef + safety increase
    xx = 0:bar_L/1000:bar_L
    strain_limit = 1.1 .* [maximum(x[1] for x in eRef.(xx,α,β));
                           minimum(x[1] for x in eRef.(xx,α,β))]
    
    dataset = create_dataset(numDataPts, x -> bar_E * x, strain_limit[2], strain_limit[1])
    
    # define the truss problem
    problem = TrussProblem(
        A,
        [0],
        connections,
        α,
        constrained_dofs,
        node_vector = node_vector,
        num_quad_pts = 2,
        force_func = x -> force_func(α, β, x, 1.0)
    )
    
    
    nonlin_result = Datasolver.directSolverNonLinearBar(problem, dataset, NR_tol = NR_tol);
    
    # nonlin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(problem, dataset, NR_tol = NR_tol);
    
    # taking results of the last ADM iter
    uh = nonlin_result.u[end]
    uxh = uh[1:2:end]
    
    eh = nonlin_result.e[end]
    sh = nonlin_result.s[end]
    
    
    # plots
    plot(xx, uRef.(xx,β), linewidth=2, linecolor=:black)
    plot!([node_vector[i][1] for i in 1:problem.num_node], uxh, linewidth=2, linecolor=:royalblue)
    
    
    plot(xx,eRef.(xx,α,β), linewidth=2, linecolor=:black)
    plot!([node_vector[i][1] for i in 1:problem.num_node], [eh[1];eh], linewidth=2,linetype=:steppre)
    
    
    plot(xx,sRef.(xx,α,β), linewidth=2, linecolor=:black)
    plot!([node_vector[i][1] for i in 1:problem.num_node], [sh[1];sh], linewidth=2,linetype=:steppre)
end

@testset "Convergence of 1D nonlinear bar with manifactured solution" begin
    bar_L = Float64(π)
    A = 2000/1e6        # [m²]
    bar_E = 1.622e+03   # [Pa] 

    α = 1.0         # 0: linear strain     1: nonlinear strain

    NR_max_iter = 50
    NR_tol = 1e-8

    β = 0.1*π

    uRef(x,β) = β * sin(π*x/bar_L)

    eRef(x,α,β) = β * π / bar_L * cos(π*x/bar_L) * ( 1 + 0.5*α * β * π / bar_L * cos(π*x/bar_L) )

    sRef(x,α,β) = bar_E * eRef(x,α,β)

    force_func(α, β, x, λ) = λ .*
        [bar_E * A * ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_L))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_L)) / bar_L))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_L)) / bar_L) + 2) / bar_L^(2)); 0]


    N_datapoints = [2^n for n in 2:9]
    N_elements = [2^n for n in 2:9]

    xx = 0:bar_L/1000:bar_L
    strain_limit = 1.1 .* [maximum(x[1] for x in eRef.(xx,α,β));
                        minimum(x[1] for x in eRef.(xx,α,β))]

    # allocation
    l2e = zeros(length(N_datapoints),length(N_elements))

    for (i,N_d) in enumerate(N_datapoints), (j,N_e) in enumerate(N_elements)
        @show N_d, N_e
        local dataset = create_dataset(N_d, x -> bar_E * x, strain_limit[2], strain_limit[1])

        # mesh
        h = bar_L/N_e
        node_vector = [ [(i-1)*h, 0] for i in 1:N_e+1 ]

        constrained_dofs = [
            (1, 1),
            (1, 2),
            (N_e+1,1),
            (N_e+1,2)
        ]

        connections = [ (i, i+1) for i in 1:N_e ]

        local nonlinear_problem = TrussProblem(
            A,
            [0],
            connections,
            α,
            constrained_dofs,
            node_vector = node_vector,
            num_quad_pts = 2,
            force_func = x -> force_func(α, β, x, 1.0)
        )

        # local result = greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = NR_tol)

        local result = directSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = NR_tol)

        l2e[i,j] = Datasolver.relL2err1D(problem=nonlinear_problem, uNodal=result.u[end], uAfunction=x->uRef(x,β))
    end

    # plots
    contour(
            N_elements, N_datapoints, log10.(l2e'),
            ylabel = "Number of data points",
            xlabel = "Number of elements",
            colorbar_title = L"Relative $L^2$ error (log10)",
            scale = :log10,
            fill = false,
            framestyle = :box,
        )
end
