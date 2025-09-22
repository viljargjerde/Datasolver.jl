using QuadGK
using Interpolations
using Plots
using Plots.PlotMeasures
using Printf


"""
	integrate(x, y) -> Float64

Integrates a set of discrete points represented by vectors `x` and `y` using a linear interpolation function.

# Arguments
- `x::Vector`: The x-coordinates of the data points.
- `y::Vector`: The y-coordinates of the data points.

# Returns
- The result of the integration of the interpolated function over the range of `x`.
"""
function integrate(x, y)
	itp = linear_interpolation(x, y)
	a, b = minimum(x), maximum(x)
	result, error = quadgk(itp, a, b)
	return result
end

"""
	create_dataset(N, strain_stress_relation, min_strain, max_strain, min_stress, max_stress; noise_magnitude=0.0) -> Dataset

Creates a synthetic dataset for strain-stress relationships with optional noise.

# Arguments
- `N::Int`: Number of data points.
- `strain_stress_relation::Function`: Function that defines the relationship between strain and stress.
- `min_strain::Float64`: Minimum strain value.
- `max_strain::Float64`: Maximum strain value.
- `min_stress::Float64`: Minimum stress value.
- `max_stress::Float64`: Maximum stress value.
- `noise_magnitude::Float64`: Magnitude of random noise to add to the stress values (default: `0.0`).

# Returns
- A `Dataset` instance containing the generated strain and stress values.
"""
function create_dataset(N, strain_stress_relation, min_strain, max_strain, min_stress, max_stress; noise_magnitude = 0.0)
	strain = collect(range(min_strain, max_strain, N))
	stress = strain_stress_relation.(strain)
	stress .-= minimum(stress) - min_stress
	stress *= (max_stress) / (maximum(stress))
	if noise_magnitude > 0
		stress .+= (randn(length(stress)) .- 0.5) .* noise_magnitude
	end
	return Dataset(strain, stress)
end

"""
	create_dataset(N, strain_stress_relation, min_strain, max_strain; noise_magnitude=0.0) -> Dataset

Creates a synthetic dataset for strain-stress relationships with optional noise.

# Arguments
- `N::Int`: Number of data points.
- `strain_stress_relation::Function`: Function that defines the relationship between strain and stress.
- `min_strain::Float64`: Minimum strain value.
- `max_strain::Float64`: Maximum strain value.
- `noise_magnitude::Float64`: Magnitude of random noise to add to the stress values (default: `0.0`).

# Returns
- A `Dataset` instance containing the generated strain and stress values.
"""
function create_dataset(N, strain_stress_relation, min_strain, max_strain; noise_magnitude = 0.0)
	strain = collect(range(min_strain, max_strain, N))
	stress = strain_stress_relation.(strain)
	if noise_magnitude > 0
		stress += randn(length(stress)) * noise_magnitude
	end
	return Dataset(strain, stress)
end

"""
	connect_in_sequence(N_elements) -> Vector{Vector{Int}}

Creates a list of connections in sequence for a set of elements.

# Arguments
- `N_elements::Int`: Number of elements to connect.

# Returns
- A vector of vectors representing sequential connections.
"""
function connect_in_sequence(N_elements)
	return [[i, i + 1] for i in 1:N_elements]
end

"""
	create_Φ_bar(N_nodes, total_length) -> Vector{Vector{Float64}}

Generates the coordinates of nodes for a 1D bar of specified length.

# Arguments
- `N_nodes::Int`: Number of nodes.
- `total_length::Float64`: The total length of the bar.

# Returns
- A vector of coordinates representing the positions of nodes along the bar.
"""
function create_Φ_bar(N_nodes, total_length)
	step_size = total_length / (N_nodes - 1)
	return [[(i - 1) * step_size, 0] for i in 1:N_nodes]
end

"""
	get_integration_interval(Φ, i) -> Tuple{Float64, Float64}

Determines the integration interval for a given node based on its neighbors.

# Arguments
- `Φ::Vector{Vector{Float64}}`: Node positions.
- `i::Int`: Index of the current node.

# Returns
- A tuple `(start, end)` representing the integration interval.
"""
function get_integration_interval(Φ, i)
	if i == 1
		before = 0.0
	else
		before = (Φ[i][1] - Φ[i-1][1]) / 2
	end
	if i == length(Φ)
		after = 0.0
	else
		after = (Φ[i+1][1] - Φ[i][1]) / 2
	end
	middle = Φ[i][1]
	return (middle - before, middle + after)
end

"""
	plot_configuration(Φ, connections) -> Plot

Plots the initial configuration of a structure defined by node positions and connections.

# Arguments
- `Φ::Vector{Vector{Float64}}`: Node positions.
- `connections::Vector{Vector{Int}}`: Connections between nodes.

# Returns
- A plot object representing the configuration.
"""
function plot_configuration(Φ, connections)
	# Extract x and y coordinates from the points in Φ
	xs = [p[1] for p in Φ]
	ys = [p[2] for p in Φ]

	plt = plot(xlabel = "x", ylabel = "y", title = "Initial configuration", legend = false)
	# Plot the connections
	for (i, connection) in enumerate(connections)
		a = Φ[connection[1]]
		b = Φ[connection[2]]
		xs_line = [a[1], b[1]]
		ys_line = [a[2], b[2]]
		if i != length(connections)
			plot!(xs_line, ys_line, label = "", color = :black, linewidth = 3)
		else
			plot!(xs_line, ys_line, label = "Elements", color = :black, linewidth = 3)
		end
	end

	# Plot nodes
	plt = scatter!(xs, ys, markersize = 6, color = :black,
		markershape = :circle, markercolor = :white, label = "Nodes")

	return plt
end

"""
	discretice_1d_force(force_function, Φ) -> Vector{Float64}

Discretizes a 1D force function by integrating over the intervals defined by node positions `Φ`.

# Arguments
- `force_function::Function`: The force function to discretize.
- `Φ::Vector{Vector{Float64}}`: Node positions.

# Returns
- A vector representing the discretized force at each node.
"""
function discretice_1d_force(force_function, Φ)
	return [quadgk(force_function, get_integration_interval(Φ, i)...)[1] for i in eachindex(Φ)]
end

"""
	setup_1d_bar_problem(N_elements, L, force_function) -> Tuple

Sets up a 1D bar problem, defining the connections, node positions, discretized forces, and fixed degrees of freedom.

# Arguments
- `N_elements::Int`: Number of elements in the bar.
- `L::Float64`: Total length of the bar.
- `force_function::Function`: The force function applied to the bar.

# Returns
- A tuple `(connections, Φ, f, fixed_dof)` representing the setup of the bar problem.
"""
function setup_1d_bar_problem(N_elements, L, force_function, remove_dofs = true)
	connections = connect_in_sequence(N_elements)
	Φ = create_Φ_bar(N_elements + 1, L)
	if remove_dofs
		f = discretice_1d_force(force_function, Φ)[begin+1:end-1] # Manually removing DOF from f
	else
		f = discretice_1d_force(force_function, Φ)
	end
	fixed_dof = [(i, 2) for i in eachindex(Φ)]
	append!(fixed_dof, [(1, 1), (N_elements + 1, 1)])
	return connections, Φ, f, fixed_dof
end

"""
	plot_dataset(dataset::Dataset; legend=false, title="Strain stress dataset") -> Plot

Plots a strain-stress dataset.

# Arguments
- `dataset::Dataset`: The dataset to plot.
- `legend::Bool`: Whether to display the legend (default: `false`).
- `title::String`: Title of the plot (default: "Strain stress dataset").

# Returns
- A plot object of the strain-stress dataset.
"""
function plot_dataset(dataset::Dataset; legend = false, title = "Strain stress dataset")
	scatter(dataset.E, dataset.S, xlabel = "Strain", ylabel = "Stress", title = title, legend = legend, markersize = 3, markercolor = :white)
end

"""
	plot_dataset(dataset::Dataset, final_results; legend=true, title="Strain stress dataset") -> Plot

Plots a strain-stress dataset alongside final computed results.

# Arguments
- `dataset::Dataset`: The original dataset.
- `final_results::SolveResults`: The final results of a solver.
- `legend::Bool`: Whether to display the legend (default: `true`).
- `title::String`: Title of the plot (default: "Strain stress dataset").

# Returns
- A plot object combining the original dataset and the solver results.
"""
function plot_dataset(dataset::Dataset, final_results; legend = true, title = "Strain stress dataset")
	scatter(dataset.E, dataset.S, xlabel = "Strain", ylabel = "Stress", title = title, label = "Dataset", markersize = 3, markercolor = :white)
	scatter!(final_results.e, final_results.s, markersize = 6, markercolor = :green, label = "Stress and strain")
	scatter!(final_results.E, final_results.S, markersize = 2.5, markercolor = :red, label = "Selected points", marker = :square, legend = legend)
end

"""
	get_rel_diff(xs, u_solved, u_analytical) -> Float64

Computes the relative difference between the solved displacement `u_solved` and analytical displacement `u_analytical`.

# Arguments
- `xs::Vector`: x-coordinates of the nodes.
- `u_solved::Vector`: Solved displacement values.
- `u_analytical::Vector`: Analytical displacement values.

# Returns
- The relative difference between `u_solved` and `u_analytical`.
"""
function get_rel_diff(xs, u_solved, u_analytical)
	return sqrt(integrate(xs, (u_solved - u_analytical) .^ 2)) / (sqrt(integrate(xs, u_analytical .^ 2)))
end

"""
	convergence_analysis(results::Vector{NamedTuple}, us)

Performs a convergence analysis by plotting the relative differences for different numbers of elements and data points.

# Arguments
- `results::Vector{NamedTuple}`: A vector of named tuples containing the results of various simulations.
- `us::Vector`: Analytical displacement values for comparison.

# Returns
- A contour plot representing the relative differences (in log10 scale).
"""
function convergence_analysis(results::Vector{NamedTuple}, us)
	all_N_d = sort!(unique([r.N_datapoints for r in results]))
	all_N_e = sort!(unique([r.N_elements for r in results]))
	rel_differences = zeros(length(all_N_e), length(all_N_d))

	for (result, u) in zip(results, us)
		xs = [p[1] for p in result.result.Φ]
		u_solved = get_final(result.result).u
		if length(u_solved) != length(xs)
			u_solved = [0.0, u_solved..., 0.0]
		end
		rel_difference = get_rel_diff(xs, u_solved, u)
		i = findfirst(all_N_e .== result.N_elements)
		j = findfirst(all_N_d .== result.N_datapoints)
		rel_differences[i, j] = rel_difference
	end

	contour(all_N_d, all_N_e, log10.(rel_differences),
		xlabel = "N datapoints", ylabel = "N elements", scale = :log2,
		title = "Relative Difference Contour Plot (Log10 Scale)", fill = true)
end

# TODO generalize to 2D
"""
	plot_results(result::SolveResults; dataset=nothing, title="") -> Plot

Plots the results of a solver, including strain, stress, displacement, and more.

# Arguments
- `result::SolveResults`: The results of a solver.
- `dataset::Union{Dataset, Nothing}`: The original dataset (default: `nothing`).
- `title::String`: Title of the plot (default: `""`).

# Returns
- A plot object representing the various solver results.
"""
function plot_results(result::SolveResults; dataset = nothing, title = "")
	# Extract x positions from result.Φ
	x_nodes = [p[1] for p in result.Φ]
	@show x_nodes

	# Calculate midpoints between nodes for lambda and mu
	x_midpoints = [(x_nodes[i] + x_nodes[i+1]) / 2 for i in 1:length(x_nodes)-1]

	# Create the subplots
	if dataset === nothing
		plot_layout = @layout [a b c; d _ e; f g h]
	else
		plot_layout = @layout [a b c; d e f; g h i]
	end
	final_result = get_final(result)
	tick_formatter = x -> @sprintf("%.2g", x)
	# Plot e, s, and u at each node

	p1 = plot(x_midpoints, final_result.e, xlabel = "x", ylabel = "e", title = "e", marker = :x, legend = false, yformatter = tick_formatter)
	p2 = plot(x_midpoints, final_result.s, xlabel = "x", ylabel = "s", title = "s", marker = :x, legend = false, yformatter = tick_formatter)
	if length(x_nodes) > length(final_result.u)
		p3 = plot(x_nodes, [0.0, final_result.u..., 0.0], xlabel = "x", ylabel = "u", title = "u", marker = :x, legend = false, yformatter = tick_formatter)
	else
		p3 = plot(x_nodes, final_result.u, xlabel = "x", ylabel = "u", title = "u", marker = :x, legend = false, yformatter = tick_formatter)
	end
	# Plot lambda and mu
	if !isempty(final_result.λ)
		if length(x_nodes) > length(final_result.λ)
			p4 = plot(x_nodes, [0.0, final_result.λ..., 0.0], xlabel = "x", ylabel = "λ", title = "λ", marker = :x, legend = false, yformatter = tick_formatter)
		else
			@show norm(final_result.λ)
			p4 = plot(x_nodes, final_result.λ, xlabel = "x", ylabel = "λ", title = "λ", marker = :x, legend = false, yformatter = tick_formatter)
		end
	else
		p4 = plot([], [], xlabel = "x", ylabel = "λ", title = "λ", legend = false)  # Empty plot if λ is empty
	end

	if !isempty(final_result.μ)
		@show norm(final_result.μ)
		p5 = plot(x_midpoints, final_result.μ, xlabel = "x", ylabel = "μ", title = "μ", marker = :x, legend = false, yformatter = tick_formatter)
	else
		p5 = plot([], [], xlabel = "x", ylabel = "μ", title = "μ", legend = false)  # Empty plot if μ is empty
	end
	# Plot cost, balance, and compatibility
	p6 = plot(1:length(result.cost), result.cost, xlabel = "Iteration", ylabel = "Cost", title = "Cost", marker = :x, legend = false, yformatter = tick_formatter)
	p7 = plot(1:length(result.balance), [sqrt(sum(result.balance[i] .^ 2)) for i in eachindex(result.balance)], xlabel = "Iteration", ylabel = "||balance||", title = "Balance eq.", marker = :x, legend = false, yformatter = tick_formatter)
	p8 = plot(
		1:length(result.compatibility),
		[sqrt(sum(result.compatibility[i] .^ 2)) for i in eachindex(result.compatibility)],
		xlabel = "Iteration",
		ylabel = "||Compatibility||",
		title = "Compatibility",
		marker = :x,
		legend = false,
		yformatter = tick_formatter,
	)
	# Combine all the plots into a grid
	if dataset === nothing
		plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = plot_layout, size = (900, 600), right_margin = 2mm, left_margin = 4mm, plot_title = title)
	else
		data_plot = plot_dataset(dataset, final_result; legend = false, title = "Dataset")
		plot(p1, p2, p3, p4, p5, data_plot, p6, p7, p8, layout = plot_layout, size = (900, 600), right_margin = 2mm, left_margin = 4mm, plot_title = title)
	end
end

