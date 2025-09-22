using LinearAlgebra
using Datasolver

pgfplotsx()
default(size = (800, 600))


# # inputs
bar_length = 1.0 * 1000      # [mm]   - initial length of the bar
area = 400.0     # [mm^2] - cross-sectional area of the bar
force = x -> [0.1]  # [kN/mm]   - constant uniform distributed load
num_ele = 6     # [-]   - number of elements
bar_E = 1e3;        # [MPa]  - Young_modulus
numDataPts = 21;    # [-]   - number of data points, odd number to ensure zero strain is in the dataset
function get_problems(num_ele)
	linear_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		0.0;
		right_fixed = true,
	)
	nonlinear_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		num_ele,
		1.0;
		right_fixed = true,
	)
	linear_problem, nonlinear_problem

end

linear_problem, nonlinear_problem = get_problems(num_ele)



# Generate data 
dist_F = linear_problem.force(nothing)
strain_limit = norm(dist_F) * linear_problem.length / (2 * bar_E * linear_problem.area);
if norm(linear_problem.force(1)) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)


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

	plot()
	line_col = names(grouped)[1]
	for line_name in unique(grouped[!, line_col])
		filtered_df = filter(row -> row[line_col] == line_name, grouped)
		x_name = names(grouped)[2]
		x = filtered_df[!, x_name]
		y = filtered_df[!, "Median $(lowercase(y_col_name))"]
		lower_err = y .- filtered_df.q25
		upper_err = filtered_df.q75 .- y

		# Cap lower error to avoid going below zero
		lower_err = map((yi, le) -> min(le, yi), y, lower_err)

		yerr = (lower_err, upper_err)
		@show yerr
		plot!(x, y, label = line_name, marker = :circle, scale = :log2,
			yerror = yerr, xlabel = x_name, ylabel = "Median $(lowercase(y_axis_name))")

	end
	savefig(replace(results_file, ".json" => ".tex"))
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

	savefig(replace(results_file, "results.json" => "plot.tex"))
end
