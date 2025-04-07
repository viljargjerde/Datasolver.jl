using LinearAlgebra
using Datasolver
using Plots.PlotMeasures
using Polynomials
using PGFPlotsX
using Plots
using Roots
using ColorSchemes
using Format
using Printf

paired_colors = colorschemes[:tableau_20]
single_colors = colorschemes[:tableau_10]
# single_colors = colorschemes[:tableau_20][2:2:end]
pgfplotsx()
default(size = (400, 300), markersize = 3, palette = single_colors)

# # inputs
bar_length = 1.0 * 1000      # [mm]   - initial length of the bar
area = 400.0     # [mm^2] - cross-sectional area of the bar
force = x -> [100.0]  # [N/mm]   - constant uniform distributed load
bar_E = 2.7e3;        # [MPa]  - Young_modulus
num_ele = 8     # [-]   - number of elements
numDataPts = 17;    # [-]   - number of data points, odd number to ensure zero strain is in the dataset

all_results_file = "../master_thesis/all_results.tex"

function update_tex_command(filename::String, cmd_name::String, value::String)
	commands = Dict{String, String}()
	other_lines = String[]  # To preserve non-command lines
	pattern = r"\\newcommand\{\\([^\}]+)\}\{(.*)\}"

	# Read file and extract existing commands
	for line in eachline(filename)
		m = match(pattern, line)
		if m !== nothing
			commands[m.captures[1]] = m.captures[2]
		else
			push!(other_lines, line)
		end
	end

	# Update or insert new command
	commands[cmd_name] = value

	# Write everything back to file
	open(filename, "w") do io
		for (cmd, val) in sort(collect(commands))
			println(io, "\\newcommand{\\$cmd}{$val}")
		end
		for line in other_lines
			println(io, line)
		end
	end
end




function get_problems(num_ele)
	linear_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		0.0;
		right_fixed = false,
	)
	nonlinear_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		1.0;
		right_fixed = false,
	)
	linear_problem, nonlinear_problem

end

linear_problem, nonlinear_problem = get_problems(num_ele)



# Generate data 
dist_F = linear_problem.force(nothing)
strain_limit = norm(dist_F) * linear_problem.length / (bar_E * linear_problem.area);
if norm(linear_problem.force(1)) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E * x, 0, 2 * strain_limit)


scatter(dataset.E, dataset.S,
	xlabel = "Strain",
	ylabel = "Stress [MPa]",
	legend = false, markercolor = nothing)

savefig("../master_thesis/figures/dataset.tex")

function check_similarity(results1, results2)
	is_same = true
	fields = [
		"E",
		"S",
		"s",
		"e",
		"u",
		"cost"]
	for field in fields
		if !(results1[field] ≈ results2[field])
			println("$field is not the same, reldiff: $(norm(results1[field] - results2[field])/norm(results1[field]))")
			is_same = false
			break
		end
	end
	is_same

end


update_tex_command(all_results_file, "barLength", string(Int(bar_length)))
update_tex_command(all_results_file, "numDataPts", string(numDataPts))
update_tex_command(all_results_file, "numEle", string(num_ele))
update_tex_command(all_results_file, "area", string(Int(area)))
update_tex_command(all_results_file, "forceConst", string(Int(force(1)[1])))
update_tex_command(all_results_file, "barE", string(Int(bar_E)))
update_tex_command(all_results_file, "minstressData", @sprintf("%.0f", minimum(dataset.S)))
update_tex_command(all_results_file, "maxstressData", @sprintf("%.0f", maximum(dataset.S)))
update_tex_command(all_results_file, "minstrainData", @sprintf("%.3f", minimum(dataset.E)))
update_tex_command(all_results_file, "maxstrainData", @sprintf("%.3f", maximum(dataset.E)))


function process_results(df, results_file, y = ("Solve time", "Solve time (s)"))
	(y_col_name, y_axis_name) = y

	cols = [n for n in names(df) if n ∉ ["Result", y_col_name]]
	# Compute statistics
	grouped = combine(groupby(df, cols),
		y_col_name => median => "Median $(lowercase(y_col_name))",
		y_col_name => (x -> quantile(x, 0.25)) => :q25,
		y_col_name => (x -> quantile(x, 0.75)) => :q75,
	)


	table = unstack(grouped, cols[2], cols[1], "Median $(lowercase(y_col_name))")

	open(replace(results_file, ".json" => ".tex"), "w") do f
		pretty_table(f, table, backend = Val(:latex), show_subheader = false)
		pretty_table(table, show_subheader = false)
	end

	# p = plot(
	# 	legend = :topleft,
	# )

	# line_col = names(grouped)[1]
	# for line_name in unique(grouped[!, line_col])
	# 	filtered_df = filter(row -> row[line_col] == line_name, grouped)
	# 	x_name = names(grouped)[2]
	# 	x = filtered_df[!, x_name]
	# 	y = filtered_df[!, "Median $(lowercase(y_col_name))"]
	# 	lower_err = y .- filtered_df.q25
	# 	upper_err = filtered_df.q75 .- y

	# 	# Cap lower error to avoid going below zero
	# 	lower_err = map((yi, le) -> min(le, yi), y, lower_err)

	# 	yerr = (lower_err, upper_err)
	# 	@show yerr
	# 	plot!(x, y, label = line_name, marker = :circle, scale = :log2,
	# 		xlabel = x_name, ylabel = "Median $(lowercase(y_axis_name))", markercolor = :transparent)

	# 	a, b, f = estimate_powerlaw(x[2:end], y[2:end])
	# 	x_fit = range(minimum(x), stop = maximum(x), length = 100)
	# 	y_fit = a .* x_fit .^ b
	# 	plot!(x_fit, y_fit, label = nothing, lw = 2, linestyle = :dash)

	# end
	# display(p)
	# savefig(replace(results_file, "results.json" => "figure.tex"))
	return table
end


function process_results_3vars(df, results_file; y = ("Solve time", "Solve time (s)"), scatter = false)
	(y_col_name, y_axis_name) = y
	cols = [n for n in names(df) if n ∉ ["Result", y_col_name]]
	third_var = cols[3]  # Third variable for splitting tables and marker shapes

	# Compute statistics
	grouped = combine(groupby(df, cols),
		y_col_name => median => "Median $(lowercase(y_col_name))",
		y_col_name => (x -> quantile(x, 0.25)) => :q25,
		y_col_name => (x -> quantile(x, 0.75)) => :q75,
	)

	# Generate tables for each unique value in the third column
	for unique_value in unique(grouped[!, third_var])
		sub_df = filter(row -> row[third_var] == unique_value, grouped)
		table = unstack(sub_df, cols[2], cols[1], "Median $(lowercase(y_col_name))")
		filename = replace(results_file, ".json" => "_$(unique_value).tex")

		open(filename, "w") do f
			pretty_table(f, table, backend = Val(:latex), show_subheader = false)
			pretty_table(table, show_subheader = false)
		end
	end

	# Define marker shapes
	markers = [:square, :diamond, :star5, :circle, :utriangle, :dtriangle]
	unique_values = unique(grouped[!, third_var])
	shape_dict = Dict(zip(unique_values, markers))  # Assign shapes

	# Plot results
	plot()
	line_col = cols[1]
	for line_name in unique(grouped[!, line_col])
		for unique_value in unique_values
			filtered_df = filter(row -> row[line_col] == line_name && row[third_var] == unique_value, grouped)
			x_name = cols[2]
			x = filtered_df[!, x_name]
			y = filtered_df[!, "Median $(lowercase(y_col_name))"]
			err = filtered_df.q75 .- filtered_df.q25
			@show x
			@show y
			plot!(x, y, label = "$(line_name) - $(unique_value)", marker = shape_dict[unique_value],
				xlabel = x_name, ylabel = "Median $y_axis_name", seriestype = scatter ? :scatter : :line)
		end
	end

	savefig(replace(results_file, "results.json" => "figure.tex"))
end

function estimate_powerlaw(x, y)
	logx = log2.(x)
	logy = log2.(y)
	p = fit(logx, logy, 1)  # logy = b * logx + log2(a)
	b = p.coeffs[2]
	a = 2^p.coeffs[1]
	println("Estimated power law: y = $a * x^$b")
	f(x) = a * x^b
	return a, b, f
end

function estimate_quadratic_powerlaw(x, y)
	logx = log2.(x)
	logy = log2.(y)
	p = fit(logx, logy, 2)  # logy = c2 * logx^2 + c1 * logx + c0
	c0, c1, c2 = p.coeffs
	println("Estimated model: log2(y) = $c2 * log2(x)^2 + $c1 * log2(x) + $c0")
	f(x) = 2^(c0 + c1 * log2(x) + c2 * log2(x)^2)
	return c0, c1, c2, f
end
