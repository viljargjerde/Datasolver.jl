using LinearAlgebra
using Interpolations
using Symbolics
using Memoization

"""
	struct Dataset{T <: Real}

Represents a dataset containing vectors of empirical strain `E`, stress `S` values, and a numerical constant `C` typically calculated using sum(S ./ (E .+ 1e-10)) / length(S).

# Fields
- `E::Vector{T}`: A vector of strain values.
- `S::Vector{T}`: A vector of stress values.
- `C::T`: A computed constant based on `E` and `S` values.
"""
struct Dataset{T<:Real}
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
    Φ::Union{Vector{Vector{Float64}},Vector{Float64}}
    e::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    E::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    s::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    S::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    u::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    λ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    μ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    cost::Vector{Float64} = Vector{Vector{Float64}}()
    balance::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
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
    return sum(S ./ (E .+ 1e-10)) / length(S)
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
	create_B_matrix(connections, Φ) -> Matrix

Creates the `B` matrix for the given set of connections and node positions `Φ`.

# Arguments
- `connections::Vector`: A vector of node connection pairs.
- `Φ::Vector{Vector{Float64}}`: Node coordinates or positions.

# Returns
- The constructed `B` matrix.
"""
function create_B_matrix(connections, Φ)
    B = zeros(length(connections), 2 * length(Φ))
    I_mat = Matrix{Float64}(I, 2, 2)

    # Iterate over connections
    for (k, (i, j)) in enumerate(connections)
        seg = Φ[j] - Φ[i]
        d_div_L = seg / norm(seg)^2
        local_B = d_div_L' * [-I_mat I_mat]

        # Mapping local to global DOFs
        glob_dofs = [2 * i - 1, 2 * i, 2 * j - 1, 2 * j]
        # Place local B into global B matrix
        B[k, glob_dofs] = local_B
    end
    return B
end

function create_B_matrix_1D(connections, Φ)
    B = zeros(length(connections), length(Φ))

    # Iterate over connections
    for (k, (i, j)) in enumerate(connections)
        seg = Φ[j] - Φ[i]
        d_div_L = 1 / seg  # In 1D, d_div_L is just 1 over the length of the element
        local_B = d_div_L * [-1.0, 1.0]

        # Mapping local to global DOFs
        glob_dofs = [i, j]

        # Place local B into global B matrix
        B[k, glob_dofs] = local_B
    end
    return B
end

"""
	remove_dofs(B::Matrix, node_coordinate_pairs) -> Matrix

Removes degrees of freedom (DOFs) from matrix `B` based on specified node-coordinate pairs.

# Arguments
- `B::Matrix`: The `B` matrix.
- `node_coordinate_pairs::Vector{Tuple}`: A vector of node-coordinate pairs, e.g., `[(2, 1), (5, 2)]`.

# Returns
- A new matrix `B` with specified DOFs removed.
"""
function remove_dofs(B::Matrix, node_coordinate_pairs)
    columns_to_remove = (2 * (n - 1) + c for (n, c) in node_coordinate_pairs)
    return B[:, setdiff(1:end, columns_to_remove)]
end

"""
	remove_dofs(f::Vector, node_coordinate_pairs) -> Vector

Removes degrees of freedom (DOFs) from vector `f` based on specified node-coordinate pairs.

# Arguments
- `f::Vector`: The input vector.
- `node_coordinate_pairs::Vector{Tuple}`: A vector of node-coordinate pairs, e.g., `[(2, 1), (5, 2)]`.

# Returns
- A new vector `f` with specified DOFs removed.
"""
function remove_dofs(f::Vector, node_coordinate_pairs)
    elems_to_remove = (2 * (n - 1) + c for (n, c) in node_coordinate_pairs)
    return f[setdiff(1:end, elems_to_remove)]
end

"""
	solve_system(B, C, w, E, S, f) -> Tuple

Solves a system of equations involving `B`, `C`, weights `w`, strain `E`, stress `S`, and force `f`.

# Arguments
- `B::Matrix`: The `B` matrix.
- `C::Vector`: Compatibility coefficients.
- `w::Vector`: Weights.
- `E::Vector`: strain values.
- `S::Vector`: Stress values.
- `f::Vector`: Force values.

# Returns
- A tuple containing `(ϵ, σ, u, λ, μ)`, which are the solution values.
"""
function solve_system(B, C, w, E, S, f)
    Z = zeros
    k, n = size(B)
    A = [
        Z(n, n) Z(n, k) Z(n, k) B' Z(n, n);
        Z(k, n) Diagonal(w .* C) Z(k, k) I Z(k, n);
        Z(k, n) Z(k, k) Diagonal(w ./ C) Z(k, k) B;
        -B I Z(k, k) Z(k, k) Z(k, n);
        Z(n, n) Z(n, k) (w .* B)' Z(n, k) Z(n, n)
    ]

    b = vcat(
        Z(n),
        w .* C .* E,
        (w ./ C) .* S,
        Z(k),
        f,
    )

    x = A \ b
    u = x[begin:n]
    ϵ = x[n+1:n+k]
    σ = x[n+k+1:n+2*k]
    μ = x[n+2*k+1:n+3*k]
    λ = x[n+3*k+1:end]

    return ϵ, σ, u, λ, μ
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
	calc_w(connections, Φ, A) -> Vector

Calculates `w` for each connection based on node positions `Φ` and cross-sectional areas `A`.

# Arguments
- `connections::Vector`: A vector of node connection pairs.
- `Φ::Vector{Vector{Float64}}`: Node coordinates or positions.
- `A::Vector`: Cross-sectional areas.

# Returns
- A vector of `w`.
"""
function calc_w(connections, Φ, A)
    return [norm(Φ[i] - Φ[j]) for (i, j) in connections] .* A
end

"""
	average_data(x::Vector, y::Vector, n::Int) -> Tuple

Averages `y` values in each of `n` intervals along `x` to produce a simplified representation of the data.

# Arguments
- `x::Vector`: Independent variable values.
- `y::Vector`: Dependent variable values.
- `n::Int`: Number of intervals for averaging.

# Returns
- A tuple `(mid_xs, mean_values)` representing the midpoints of intervals and averaged values.
"""
function average_data(x::Vector, y::Vector, n::Int)
    # Ensure x and y have the same length
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end

    # Get the range of x
    x_min = minimum(x)
    x_max = maximum(x)

    # Define the interval size
    interval_size = (x_max - x_min) / n

    # Initialize a vector to store the mean of y in each interval
    new_ys = Vector{Float64}()
    new_xs = Vector{Float64}()

    for i in 1:n
        # Define the boundaries of the current interval
        x_lower = x_min + (i - 1) * interval_size
        x_upper = x_min + i * interval_size
        # Get indices of x values that fall within the current interval
        indices = (x_lower .<= x) .& (x .< x_upper)

        # Calculate the mean of x and y values for the current interval
        if any(indices) # Ensure there are elements in the interval
            push!(new_xs, sum(x[indices]) / length(x[indices]))
            push!(new_ys, sum(y[indices]) / length(y[indices]))
        end
    end

    return new_xs, new_ys
end

"""
	my_new_solver(connections, Φ, A, data, f, fixed_dofs; n_smoothing_intervals=10000, verbose=false) -> SolveResults

Solves the system using given parameters and data, returning results in a `SolveResults` structure.

# Arguments
- `connections::Vector`: Connections between nodes.
- `Φ::Vector`: Node coordinates.
- `A::Vector`: Cross-sectional areas.
- `data::Dataset`: Input dataset.
- `f::Vector`: Force values.
- `fixed_dofs::Vector`: Fixed degrees of freedom.
- `n_smoothing_intervals::Int`: Number of intervals for averaging (default: `10000`).
- `verbose::Bool`: If true, outputs additional information.

# Returns
- A `SolveResults` instance containing the solution.
"""
function my_new_solver(connections, Φ, A, data, f, fixed_dofs; n_smoothing_intervals=10000, verbose=false)
    results = SolveResults(N_datapoints=length(data), Φ=Φ)
    B = create_B_matrix(connections, Φ)
    B = remove_dofs(B, fixed_dofs)
    if length(f) != size(B, 2)
        f = remove_dofs(f, fixed_dofs)
    end
    data_func = LinearInterpolation(average_data(data.S, data.E, n_smoothing_intervals)...)

    w = calc_w(connections, Φ, A)

    s = (B' \ f) ./ w
    e = data_func.(s)
    u = B \ e

    balance = B' * (s .* w) - f
    compatibility = e - B * u

    push!(results.e, e)
    push!(results.s, s)
    push!(results.E, e)
    push!(results.S, s)
    push!(results.u, u)
    push!(results.λ, zero(u))
    push!(results.μ, zero(s))
    push!(results.balance, balance)
    push!(results.compatibility, compatibility)
    push!(results.cost, 0.0)
    return results
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

"""
	datasolve(connections, Φ, A, data::Dataset, f, fixed_dofs; initialization=true, max_iterations=1000, verbose=true) -> Union{SolveResults, Nothing}

Solves the dataset optimization problem using given parameters and returns a `SolveResults` instance.

# Arguments
- `connections::Vector`: Connections between nodes.
- `Φ::Vector`: Node coordinates.
- `A::Vector`: Cross-sectional areas.
- `data::Dataset`: Input dataset.
- `f::Vector`: Force values.
- `fixed_dofs::Vector`: Fixed degrees of freedom.
- `initialization::Bool`: Whether to initialize strain and stress values (default: `true`).
- `max_iterations::Int`: Maximum number of iterations allowed (default: `1000`).
- `verbose::Bool`: If true, outputs additional information during solving.

# Returns
- A `SolveResults` instance if successful, otherwise `nothing`.
"""
function datasolve(connections, Φ, A, data::Dataset, f, fixed_dofs; initialization=true, max_iterations=1000, verbose=true)
    ########## TODO consider moving into a new "model" type: #############
    B = create_B_matrix(connections, Φ)
    n_dofs_to_fix = size(nullspace(B), 2)
    @assert n_dofs_to_fix <= length(fixed_dofs) "Problem has $n_dofs_to_fix degrees of freedom that needs fixing, but $(length(fixed_dofs)) fixed degrees of freedom were provided"
    B = remove_dofs(B, fixed_dofs)
    if length(f) != size(B, 2)
        f = remove_dofs(f, fixed_dofs)
    end
    results = SolveResults(N_datapoints=length(data), Φ=Φ)
    m = size(B, 1)
    w = calc_w(connections, Φ, A)
    ##############################################################################

    # Capital letters refer to variables from the dataset

    # i) Initialization:
    if initialization
        s = (B' \ f) ./ w
        S = zeros(size(s))
        E = zeros(size(s))
        for i in 1:m
            # choose E, S pair where S is closest to s
            best_idx = argmin((abs(data[j][2] - s[i]) for j in 1:length(data)))
            E[i] = data.E[best_idx]
            S[i] = data.S[best_idx]
        end
    else
        # Choose random E, S pair 
        points = (data[rand(1:length(data))] for _ in 1:m)
        E, S = unpack_pairs(points)
    end

    # Keep track of the history of E and S
    push!(results.E, E)
    push!(results.S, S)
    for k in 1:max_iterations

        # ii + iii) 
        e, s, u, λ, μ = solve_system(B, data.C, w, E, S, f)
        e = B * u

        # iv) 
        (E, S), cost = choose_closest_to(e, s, w, data)

        balance = B' * (s .* w) - f
        compatibility = e - B * u

        push!(results.e, e)
        push!(results.s, s)
        push!(results.E, E)
        push!(results.S, S)
        push!(results.u, u)
        push!(results.λ, λ)
        push!(results.μ, μ)
        push!(results.balance, balance)
        push!(results.compatibility, compatibility)
        push!(results.cost, cost)

        # v) 
        if results.E[end-1] == results.E[end] && results.S[end-1] == results.S[end]
            if verbose
                println("Converged in $k iterations")
            end
            if length(data) > 4 && !check_dataset_is_safe(results.E[end], results.S[end], data)
                println("WARNING dataset might be too small")
            end
            return results
        end
    end
    println("Failed to converge in $max_iterations iterations")
    return nothing
end

include("LP_solver.jl")

##################### Non linear #################################


function newton_raphson(x, g, J, free_dofs; tol=1e-6, max_iter=500)

    # Combine fixed degrees of freedom
    g_vec = g(x)
    @show norm(g_vec)
    for i in 1:max_iter
        # Solve for the update step Δx
        J_mat = J(x)
        # Remove rows and columns corresponding to fixed degrees of freedom
        @show cond(J_mat[free_dofs, free_dofs])
        Δx = -J_mat[free_dofs, free_dofs] \ g_vec[free_dofs]

        # Update only the free degrees of freedom
        # x[free_dofs] = x[free_dofs] + 0.1Δx
        x[free_dofs] = x[free_dofs] + Δx

        g_vec = g(x)

        @show norm(Δx), Δx
        @show x[free_dofs]
        @show norm(g_vec), g_vec

        # Check for convergence
        if norm(Δx) < tol && norm(g(x)) < tol
            println("Converged in $i iterations.")
            return x
        end
    end
    println("Newton-Raphson did not converge within $max_iter iterations.")
    return x
end

function get_free_dofs(fixed_dofs, x_sizes)
    fixed_dofs = vcat(fixed_dofs, fixed_dofs .+ sum(x_sizes[begin:end-1]))
    return setdiff(1:sum(x_sizes), fixed_dofs)


end


function construct_N(ϕ)

    n = length(ϕ)  # Number of nodes
    N = Vector{Num}(undef, n)
    @variables x
    for i in 1:n
        b = ϕ[i]
        if i == 1
            c = ϕ[i+1]
            slope = -1 / (c - b)
            N[i] = ifelse((x >= b) & (x <= c), 1 + slope * (x - b), 0.0)
        elseif i == n
            a = ϕ[i-1]
            # Last node: derivative is constant over [ϕ[n-1], ϕ[n]]
            slope = 1 / (b - a)
            N[i] = ifelse((x >= a) & (x <= b), 1 + slope * (x - b), 0.0)
        else
            a = ϕ[i-1]
            c = ϕ[i+1]
            # Interior nodes: piecewise constant derivative
            slope_left = 1 / (b - a)      # Slope over [a, b]
            slope_right = -1 / (c - b)    # Slope over [b, c]
            N[i] = ifelse((x >= a) & (x < b), 1 + slope_left * (x - b), 0.0) +
                   ifelse((x >= b) & (x <= c), 1 + slope_right * (x - b), 0.0)
        end
    end
    return N'
end
function construct_Np(ϕ)

    n = length(ϕ)  # Number of nodes
    Np = Vector{Num}(undef, n)
    @variables x
    for i in 1:n
        b = ϕ[i]
        if i == 1
            c = ϕ[i+1]
            # First node: derivative is constant over [ϕ[1], ϕ[2]]
            slope = -1 / (c - b)
            Np[i] = ifelse((x >= b) & (x <= c), slope, 0.0)
        elseif i == n
            a = ϕ[i-1]
            # Last node: derivative is constant over [ϕ[n-1], ϕ[n]]
            slope = 1 / (b - a)
            Np[i] = ifelse((x >= a) & (x <= b), slope, 0.0)
        else
            a = ϕ[i-1]
            c = ϕ[i+1]
            # Interior nodes: piecewise constant derivative
            slope_left = 1 / (b - a)      # Slope over [a, b]
            slope_right = -1 / (c - b)    # Slope over [b, c]
            Np[i] = ifelse((x >= a) & (x < b), slope_left, 0.0) +
                    ifelse((x >= b) & (x <= c), slope_right, 0.0)
        end
    end
    return Np'
end

# Function to construct the vector of element-based constant functions R
function construct_R(Φ, connections)

    @variables x
    n = length(connections)
    R = Vector{Num}(undef, n)

    for (i, j) in connections
        a = min(Φ[i], Φ[j])
        b = max(Φ[i], Φ[j])
        # Each function is constant (1) over the element between phi[i] and phi[j]
        R[i] = ifelse((x >= a) & (x <= b), 1.0, 0.0)
    end
    return R'
end



# General function that creates a function from a symbolic expression and integrates it from a to b
@memoize function integrate(integrand::Num, a, b)
    if iszero(integrand == 0)
        return 0.0
    end
    @variables x
    return quadgk(build_function(integrand, x, expression=Val{false}), a, b)[1]

end

# The following functions integrate all elements and adds it to the corresponding element in the matrix.
# For example int_NR_ corresponds to \int Np'*R*(integrand), where integrand is a scalar, for example R*s
# These functions are what the int_NR etc. functions inside g() and J() are calling

@memoize function int_NR_(Np, R, connections, Φ, integrand)
    result = zeros(length(Np), length(R))
    for (R_i, (i, j)) in enumerate(connections)
        a = Φ[i]
        b = Φ[j]
        result[i, R_i] += integrate(integrand * R[R_i] * Np[i], a, b)
        result[j, R_i] += integrate(integrand * R[R_i] * Np[j], a, b)
    end
    return result
end

@memoize function int_NN_(Np, connections, Φ, integrand)
    result = zeros(length(Np), length(Np))
    for (i, j) in connections
        a = Φ[i]
        b = Φ[j]
        result[i, j] += integrate(integrand * Np[i] * Np[j], a, b)
        result[j, i] += integrate(integrand * Np[i] * Np[j], a, b)
    end
    return result
end
@memoize function int_N_(N, connections, Φ, integrand)
    result = zeros(length(N))
    for (i, j) in connections
        a = Φ[i]
        b = Φ[j]
        result[i] += integrate(integrand * N[i], a, b)
        result[j] += integrate(integrand * N[j], a, b)
    end
    return result
end

@memoize function int_RR_(R, connections, Φ, integrand)
    result = zeros(length(R), length(R))
    for (R_i, (i, j)) in enumerate(connections)
        a = Φ[i]
        b = Φ[j]
        d = norm(a - b)

        # result[R_i, R_i] += integrate(integrand * R[R_i] * R[R_i], a, b) / d
        result[R_i, R_i] += (b - a) * integrand / d # Equivalent and faster
    end
    return result
end

@memoize function int_RN_(R, Np, connections, Φ, integrand)
    result = zeros(length(R), length(Np))
    for (R_i, (i, j)) in enumerate(connections)
        a = Φ[i]
        b = Φ[j]
        result[R_i, i] += integrate(integrand * R[R_i] * Np[i], a, b)
        result[R_i, j] += integrate(integrand * R[R_i] * Np[j], a, b)
    end
    return result

end
# Simple multiplication, as the integrand is not dependent on the integral variable in this case
function int_L_(connections, Φ, integrand)
    result = zeros(length(connections), length(connections))
    for (R_i, (i, j)) in enumerate(connections)
        a = Φ[i]
        b = Φ[j]
        result[R_i, R_i] += integrand * (b - a)
    end
    return result
end

function A(u, Np, R, s, C, area, connections, Φ)
    #! Changes: Only mulitplied int_L and row 5 col 3 by A 
    #! Added - to int_RR and divided by elementwise L (division inside function)
    #! Added - to row 4 col 1
    #! divided row 1 col 4, row 3 col 5 and row 4 col 1 with elementwise L

    m = length(Np)
    n = length(R)
    Z = zeros

    Ls = [norm(Φ[i] - Φ[j]) for (i, j) in connections] # Length of each element

    int_NN(integrand) = int_NN_(Np, connections, Φ, integrand)
    int_NR(integrand) = int_NR_(Np, R, connections, Φ, integrand)

    int_L(integrand) = area * int_L_(connections, Φ, integrand)
    int_RR(integrand) = -int_RR_(R, connections, Φ, integrand)
    int_RN(integrand) = int_RN_(R, Np, connections, Φ, integrand)

    # A_mat = [
    # Z(m,m)	Z(m,n)	Z(m,n)	(m,n) 	(m,m)
    # Z(n,m)	I(n,n)	Z(n,n)	(n,n)	Z(n,m)
    # Z(n,m)	Z(n,n)	I(n,n)	Z(n,n)	(n,m)
    # (n,m)		(n,n)	Z(n,n)	Z(n,n)	Z(n,m)
    # Z(m,m)	Z(m,n)	(m,n)	Z(m,n)	Z(m,m)	
    # ]
    return [
        Z(m, m) Z(m, n) Z(m, n) int_NR(1 + Np * u)./Ls' int_NN(R * s)
        Z(n, m) int_L(C) Z(n, n) -int_RR(1) Z(n, m)
        Z(n, m) Z(n, n) int_L(1 / C) Z(n, n) int_RN(1 + Np * u)./Ls
        -int_RN(1 + 1 / 2 * u' * Np')./Ls -int_RR(1) Z(n, n) Z(n, n) Z(n, m)
        Z(m, m) Z(m, n) area*int_NR((1 + Np * u)') Z(m, n) Z(m, m)
    ]

end


function g(x, N, Np, R, C, E, S, f, area, connections, Φ)
    m = length(Np)
    n = length(R)
    u, e, s, μ, λ = unpack_vector(x, [m, n, n, n, m])

    # b = [
    # Z(m)
    # n
    # n
    # Z(n)
    # m
    # ]
    Z = zeros
    int_N(integrand) = 2 * int_N_(N, connections, Φ, integrand)
    int_L(integrand) = area * int_L_(connections, Φ, integrand)
    A_mat = A(u, Np, R, s, C, area, connections, Φ)

    b = [Z(m)
        int_L(C) * E
        int_L(1 / C) * S
        Z(n)
        int_N(f)]
    res = A_mat * x - b

    return res
end

function extract_matrix_block(mat, row_sizes, col_sizes, block_row, block_col)

    end_row_idx = sum(row_sizes[begin:block_row])
    start_row_idx = end_row_idx - row_sizes[block_row] + 1
    end_col_idx = sum(col_sizes[begin:block_col])
    start_col_idx = end_col_idx - col_sizes[block_col] + 1
    return mat[start_row_idx:end_row_idx, start_col_idx:end_col_idx]
end

function J(x, Np, R, C, A, connections, Φ)
    #! Changes: Only mulitplied int_L and row 5 col 3 by A 
    #! Added - to int_RR and divided by elementwise L (division inside function)
    #! Added - to row 4 col 1
    #! divided row 1 col 4, row 3 col 5 and row 4 col 1 with elementwise L


    m = length(Np)
    n = length(R)
    u, e, s, μ, λ = unpack_vector(x, [m, n, n, n, m])
    Ls = [norm(Φ[i] - Φ[j]) for (i, j) in connections]
    Z = zeros

    int_NN(integrand) = int_NN_(Np, connections, Φ, integrand)
    int_NR(integrand) = int_NR_(Np, R, connections, Φ, integrand)
    int_L(integrand) = A * int_L_(connections, Φ, integrand)
    int_RR(integrand) = -int_RR_(R, connections, Φ, integrand)
    int_RN(integrand) = int_RN_(R, Np, connections, Φ, integrand)
    # [
    # 	(m,m)	 	Z(m, n)	(m,n)		(m,n)		(m,m)
    # 	Z(n, m) 	I(n, n)	Z(n, n)		(n,n)		Z(n,m)
    # 	(n,m)		Z(n, n)	I(n,n)		Z(n, n) 	(n,m)
    # 	(n,m) 		(n,n) 	Z(n,n)		Z(n, n)		Z(n,m)
    # 	(m,m)	 	Z(m,n) 	(m,n)	 	Z(m,n) 		Z(m,m)
    # ]

    res =
        [
            int_NN(R * μ) Z(m, n) int_NR(Np * λ) int_NR(1 + Np * u)./Ls' int_NN((R * s))
            Z(n, m) int_L(C) Z(n, n) -int_RR(1) Z(n, m)
            int_RN(Np * λ) Z(n, n) int_L(1 / C) Z(n, n) int_RN(1 + (Np * u)')./Ls
            -int_RN(1 + u' * Np')./Ls -int_RR(1) Z(n, n) Z(n, n) Z(n, m)
            int_NN(R * s) Z(m, n) A*int_NR(1 + (Np * u)') Z(m, n) Z(m, m)
        ]

    return res

end

function J_fd(x, Np, R, C, area, connections, Φ, ϵ=1e-5)

    m = length(Np)
    n = length(R)

    u, e, s, μ, λ = unpack_vector(x, [m, n, n, n, m])
    A_mat = A(u, Np, R, s, C, area, connections, Φ)
    A_mat_eps = A(u .+ ϵ, Np, R, s .+ ϵ, C, area, connections, Φ)
    return A_mat + (A_mat_eps - A_mat) ./ ϵ


end



function unpack_vector(vector, sizes)
    result = Vector{typeof(vector)}(undef, length(sizes))
    vec_i = 1
    for (i, s) in enumerate(sizes)
        result[i] = vector[vec_i:vec_i+s-1]
        vec_i += s
    end
    return result

end


function nonlin_datasolve(connections, Φ, A, data::Dataset, f, L, fixed_dofs; initialization=true, max_iterations=1000, verbose=true)
    """
    Steps: 
    Initialize, either randomly or using the same approach as earlier. 
    	Can I use the same initailization with B, or will there be a formulation using R and N?

    Solve system of equations using Newton Raphson to obtain e,s,u λ and μ
    Find closest data points in the dataset
    Check for convergence, and repeat.

    """
    B = create_B_matrix_1D(connections, Φ)
    # B = remove_dofs_1D(B, fixed_dofs)
    Np = construct_Np(Φ)#[2:end]'
    N = construct_N(Φ)
    R = construct_R(Φ, connections)
    m = length(Np) # hat / nodes 
    n = length(R) # bar / elements
    x_sizes = [m, n, n, n, m]
    free_dofs = get_free_dofs(fixed_dofs, x_sizes)

    results = SolveResults(N_datapoints=length(data), Φ=Φ)
    ### Initiialization ###
    # Choose random E, S pair 
    if initialization

        # s = (integrate.(B' * R * R', 0, L) \ int_N_(N, connections, Φ, f))
        s = integrate.(B' .* R, 0, L) \ int_N_(N, connections, Φ, f)
        S = zeros(size(s))
        E = zeros(size(s))
        for i in 1:n
            # choose E, S pair where S is closest to s
            @show best_idx = argmin((abs(data[j][2] - s[i]) for j in 1:length(data)))
            E[i] = data.E[best_idx]
            S[i] = data.S[best_idx]
        end
    else
        # Choose random E, S pair 
        points = (data[rand(1:length(data))] for _ in 1:n)
        E, S = unpack_pairs(points)
    end

    # Keep track of the history of E and S
    push!(results.E, E)
    push!(results.S, S)



    x = vcat([zeros(s) for s in x_sizes]...)
    # @show u, e, s, μ, λ = unpack_vector(x, x_sizes)

    # Create wrapper functions for more convenient use inside Newton Raphson
    g_x(x) = g(x, N, Np, R, data.C, E, S, f, A, connections, Φ)

    J_x(x) = J_fd(x, Np, R, data.C, A, connections, Φ, 1e-10)
    # J_x(x) = J(x, Np, R, data.C, A, connections, Φ)

    for i in 1:max_iterations

        # ii + iii) 
        ### Solve system of equation using Newton Raphson ###
        x = newton_raphson(x, g_x, J_x, free_dofs, max_iter=50)
        u, e, s, μ, λ = unpack_vector(x, [m, n, n, n, m])


        # iv)
        ### Choose closest match from the data ### 
        # TODO choose_closest_to needs a new version that uses the correct balance and compatibility eqs
        (E, S), cost = choose_closest_to(e, s, w, data)

        # TODO these need to be integrated before we can record them
        # balance = B' * R * s - N' * f # ? * A?
        # compatibility = e - B * u # ? R?

        push!(results.e, e)
        push!(results.s, s)
        push!(results.E, E)
        push!(results.S, S)
        push!(results.u, u)
        push!(results.λ, λ)
        push!(results.μ, μ)
        # push!(results.balance, balance)
        # push!(results.compatibility, compatibility)
        push!(results.cost, cost)

        # v) 
        if results.E[end-1] == results.E[end] && results.S[end-1] == results.S[end]
            if verbose
                println("Converged in $i iterations")
            end
            if length(data) > 4 && !check_dataset_is_safe(results.E[end], results.S[end], data)
                println("WARNING dataset might be too small")
            end
            return results
        end
    end
    println("Failed to converge in $max_iterations iterations")
    return nothing
end

