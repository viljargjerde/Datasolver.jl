using LinearAlgebra
using Interpolations



"""
	struct Dataset{T <: Real}

Represents a dataset containing vectors of empirical strain `E`, stress `S` values, and a numerical constant `C` typically calculated using sum(S ./ (E .+ 1e-10)) / length(S).

# Fields
- `E::Vector{T}`: A vector of strain values.
- `S::Vector{T}`: A vector of stress values.
- `C::T`: A computed constant based on `E` and `S` values.
"""
struct Dataset{T <: Real}
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
	Φ::Vector{Vector{Float64}}
	e::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	E::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	s::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	S::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	u::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	η::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	μ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	cost::Vector{Float64} = Vector{Vector{Float64}}()
	balance::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	compatability::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
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
	values = [getfield(result, field)[end] for field in names]
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
- A tuple containing `(ϵ, σ, u, η, μ)`, which are the solution values.
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
	η = x[n+3*k+1:end]

	return ϵ, σ, u, η, μ
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
function my_new_solver(connections, Φ, A, data, f, fixed_dofs; n_smoothing_intervals = 10000, verbose = false)
	results = SolveResults(N_datapoints = length(data), Φ = Φ)
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
	compatability = e - B * u

	push!(results.e, e)
	push!(results.s, s)
	push!(results.E, e)
	push!(results.S, s)
	push!(results.u, u)
	push!(results.η, zero(u))
	push!(results.μ, zero(s))
	push!(results.balance, balance)
	push!(results.compatability, compatability)
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
function datasolve(connections, Φ, A, data::Dataset, f, fixed_dofs; initialization = true, max_iterations = 1000, verbose = true)
	B = create_B_matrix(connections, Φ)
	B = remove_dofs(B, fixed_dofs)
	if length(f) != size(B, 2)
		f = remove_dofs(f, fixed_dofs)
	end
	results = SolveResults(N_datapoints = length(data), Φ = Φ)
	m = size(B, 1)
	w = calc_w(connections, Φ, A)

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
		e, s, u, η, μ = solve_system(B, data.C, w, E, S, f)
		e = B * u

		# iv) 
		(E, S), cost = choose_closest_to(e, s, w, data)

		balance = B' * (s .* w) - f
		compatability = e - B * u

		push!(results.e, e)
		push!(results.s, s)
		push!(results.E, E)
		push!(results.S, S)
		push!(results.u, u)
		push!(results.η, η)
		push!(results.μ, μ)
		push!(results.balance, balance)
		push!(results.compatability, compatability)
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

