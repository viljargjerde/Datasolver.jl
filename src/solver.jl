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
	Φ::Union{Vector{Vector{Float64}}, Vector{Float64}}
	e::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	E::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	s::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	S::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	u::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	λ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	μ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
	cost::Vector{Float64} = Vector{Vector{Float64}}()
	equilibrium::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
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
	filter = abs.(E) .> 1e-7
	return sum(S[filter] ./ E[filter]) / sum(filter)
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


# TODO compare to assignLocalState and choose the best implementation
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


include("LP_solver.jl")
