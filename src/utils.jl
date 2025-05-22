using Plots
using Plots.PlotMeasures
using Printf
using LaTeXStrings


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
		σ0 = noise_magnitude / 2 * maximum(stress)
		stress .+= randn(length(stress)) .* σ0                           # additive noise
		stress .+= randn(length(stress)) .* (stress .* noise_magnitude)  # proportional noise
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
		σ0 = noise_magnitude / 2 * maximum(stress)
		stress .+= randn(length(stress)) .* σ0                           # additive noise
		stress .+= randn(length(stress)) .* (stress .* noise_magnitude)  # proportional noise
	end
	return Dataset(strain, stress)
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



using Dierckx

function integrate(x, y)
	if length(x) <= 3  # too few points for cubic spline, falling back to trapezoidal rule
		return sum(diff(x) .* (y[1:end-1] .+ y[2:end]) ./ 2)
	else
		spl = Spline1D(x, y)
		return Dierckx.integrate(spl, x[1], x[end])
	end
end



"""
	calc_reldiff(xs, u_solved, u_analytical) -> Float64

Computes the relative difference between the solved displacement `u_solved` and analytical displacement `u_analytical`.

# Arguments
- `xs::Vector`: x-coordinates of the nodes.
- `u_solved::Vector`: Solved displacement values.
- `u_analytical::Vector`: Analytical displacement values.

# Returns
- The relative difference between `u_solved` and `u_analytical`.
"""
function calc_reldiff(xs, u_solved, u_analytical)
	return sqrt(integrate(xs, (u_solved - u_analytical) .^ 2)) / (sqrt(integrate(xs, u_analytical .^ 2)))
end

"""
	convergence_analysis(results::Vector{NamedTuple}, us)

Performs a convergence analysis by plotting the relative differences for different numbers of elements and data points.

# Arguments
- `results::Vector{NamedTuple}`: A vector of named tuples containing the results of various simulations.
- `us::Vector`: Analytical displacement values for comparison.

# Returns
- A contour plot representing the relative differences (in log2 scale).
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
		rel_difference = calc_reldiff(xs, u_solved, u)
		i = findfirst(all_N_e .== result.N_elements)
		j = findfirst(all_N_d .== result.N_datapoints)
		rel_differences[i, j] = rel_difference
	end
	@show rel_differences
	# contour(all_N_d, all_N_e, rel_differences,
	# 	xlabel = "N datapoints", ylabel = "N elements", scale = :log2,
	# 	title = "Relative Difference Contour Plot", fill = false)
	contour(
		all_N_d, all_N_e, log2.(rel_differences),
		xlabel = "Number of data points",
		ylabel = "Number of elements",
		colorbar_title = L"$log_2$(Relative error)",
		scale = :log2,
		fill = false,
		framestyle = :box,
	)

end

# function convergence_analysis_dual(results1::Vector{NamedTuple}, us1,
# 	results2::Vector{NamedTuple}, us2)

# 	# Extract and sort unique parameters
# 	all_N_d = sort!(unique([r.N_datapoints for r in results1]))
# 	all_N_e = sort!(unique([r.N_elements for r in results1]))

# 	rel_diff_1 = zeros(length(all_N_e), length(all_N_d))
# 	rel_diff_2 = zeros(length(all_N_e), length(all_N_d))

# 	# Fill matrix for dataset 1
# 	for (result, u) in zip(results1, us1)
# 		xs = [p[1] for p in result.result.Φ]
# 		u_solved = get_final(result.result).u
# 		if length(u_solved) != length(xs)
# 			u_solved = [0.0, u_solved..., 0.0]
# 		end
# 		i = findfirst(all_N_e .== result.N_elements)
# 		j = findfirst(all_N_d .== result.N_datapoints)
# 		rel_diff_1[i, j] = calc_reldiff(xs, u_solved, u)
# 	end

# 	# Fill matrix for dataset 2
# 	for (result, u) in zip(results2, us2)
# 		xs = [p[1] for p in result.result.Φ]
# 		u_solved = get_final(result.result).u
# 		if length(u_solved) != length(xs)
# 			u_solved = [0.0, u_solved..., 0.0]
# 		end
# 		i = findfirst(all_N_e .== result.N_elements)
# 		j = findfirst(all_N_d .== result.N_datapoints)
# 		rel_diff_2[i, j] = calc_reldiff(xs, u_solved, u)
# 	end

# 	# Shared color scale
# 	all_data = log2.([rel_diff_1; rel_diff_2])
# 	clim = extrema(all_data)

# 	plt1 = contour(
# 		all_N_d, all_N_e, log2.(rel_diff_1),
# 		xlabel = "Number of data points", ylabel = "Number of elements",
# 		clim = clim, fill = false, framestyle = :box, colorbar = false,
# 		scale = :log2,
# 		# title = "Linear strain measure",
# 	)

# 	plt2 = contour(
# 		all_N_d, all_N_e, log2.(rel_diff_2),
# 		xlabel = "Number of data points",
# 		clim = clim, fill = false, framestyle = :box,
# 		colorbar_title = L"\log_{2}(\mathrm{Relative\ error})",
# 		scale = :log2,
# 		# title = "Nonlinear strain measure",
# 	)

# 	@show minimum(rel_diff_1), minimum(rel_diff_2)
# 	plot(plt1, plt2, layout = (1, 2), link = :all)
# end
function convergence_analysis_dual(results1::Vector{NamedTuple}, us1,
	results2::Vector{NamedTuple}, us2)

	# Extract and sort unique parameters
	all_N_d = sort!(unique([r.N_datapoints for r in results1]))
	all_N_e = sort!(unique([r.N_elements for r in results1]))

	rel_diff_1 = zeros(length(all_N_e), length(all_N_d))
	rel_diff_2 = zeros(length(all_N_e), length(all_N_d))

	# Fill matrix for dataset 1
	for (result, u) in zip(results1, us1)
		xs       = [p[1] for p in result.result.Φ]
		u_solved = get_final(result.result).u
		if length(u_solved) != length(xs)
			u_solved = [0.0, u_solved..., 0.0]
		end
		i = findfirst(all_N_e .== result.N_elements)
		j = findfirst(all_N_d .== result.N_datapoints)
		rel_diff_1[i, j] = calc_reldiff(xs, u_solved, u)
	end

	# Fill matrix for dataset 2
	for (result, u) in zip(results2, us2)
		xs       = [p[1] for p in result.result.Φ]
		u_solved = get_final(result.result).u
		if length(u_solved) != length(xs)
			u_solved = [0.0, u_solved..., 0.0]
		end
		i = findfirst(all_N_e .== result.N_elements)
		j = findfirst(all_N_d .== result.N_datapoints)
		rel_diff_2[i, j] = calc_reldiff(xs, u_solved, u)
	end

	# Shared color scale
	all_data = log2.([rel_diff_1; rel_diff_2])
	clim     = extrema(all_data)

	plt1 = contour(
		all_N_d, all_N_e, log2.(rel_diff_1),
		xlabel     = "Number of data points",
		ylabel     = "Number of elements",
		clim       = clim,
		fill       = false,
		framestyle = :box,
		colorbar   = false,
		scale      = :log2,
		# title  = "Linear strain measure",
	)

	plt2 = contour(
		all_N_d, all_N_e, log2.(rel_diff_2),
		xlabel         = "Number of data points",
		clim           = clim,
		fill           = false,
		framestyle     = :box,
		colorbar_title = L"\log_{2}(\mathrm{Relative\ error})",
		scale          = :log2,
		# title        = "Nonlinear strain measure",
	)

	@show minimum(rel_diff_1), minimum(rel_diff_2)

	# <-- layout is now 2 rows × 1 column
	plot(plt1, plt2, layout = (2, 1), link = :all)
end

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
	if length(x_nodes) > length(final_result.u) # We have removed dofs
		p3 = plot(x_nodes, [0.0, final_result.u..., 0.0], xlabel = "x", ylabel = "u", title = "u", marker = :x, legend = false, yformatter = tick_formatter)
	elseif length(x_nodes) < length(final_result.u) # They are not the same dimension
		p3 = plot(x_nodes, [u[1] for u in final_result.u], xlabel = "x", ylabel = "u", title = "u", marker = :x, legend = false, yformatter = tick_formatter)
	else
		p3 = plot(x_nodes, final_result.u, xlabel = "x", ylabel = "u", title = "u", marker = :x, legend = false, yformatter = tick_formatter)
	end
	# Plot lambda and mu
	if !isempty(final_result.λ)
		if length(x_nodes) > length(final_result.λ)
			p4 = plot(x_nodes, [0.0, final_result.λ..., 0.0], xlabel = "x", ylabel = "λ", title = "λ", marker = :x, legend = false, yformatter = tick_formatter)
		elseif length(x_nodes) < length(final_result.λ)
			p4 = plot(x_nodes, [λ[1] for λ in final_result.λ], xlabel = "x", ylabel = "λ", title = "λ", marker = :x, legend = false, yformatter = tick_formatter)
		else
			p4 = plot(x_nodes, final_result.λ, xlabel = "x", ylabel = "λ", title = "λ", marker = :x, legend = false, yformatter = tick_formatter)
		end
	else
		p4 = plot([], [], xlabel = "x", ylabel = "λ", title = "λ", legend = false)  # Empty plot if λ is empty
	end

	if !isempty(final_result.μ)
		p5 = plot(x_midpoints, final_result.μ, xlabel = "x", ylabel = "μ", title = "μ", marker = :x, legend = false, yformatter = tick_formatter)
	else
		p5 = plot([], [], xlabel = "x", ylabel = "μ", title = "μ", legend = false)  # Empty plot if μ is empty
	end
	# Plot cost, equilibrium, and compatibility
	p6 = plot(1:length(result.cost), result.cost, xlabel = "Iteration", ylabel = "Cost", title = "Cost", marker = :x, legend = false, yformatter = tick_formatter)
	p7 = plot(1:length(result.equilibrium), [norm(equilibrium) for equilibrium in result.equilibrium], xlabel = "Iteration", ylabel = "||equilibrium||", title = "Equilibrium eq.", marker = :x, legend = false, yformatter = tick_formatter)
	p8 = plot(
		1:length(result.compatibility),
		[norm(c) for c in result.compatibility],
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

