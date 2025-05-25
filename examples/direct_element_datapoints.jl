include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using StatsBase: Histogram, fit



results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")

max_elements = 33
max_datapoints = 256

if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	results_list = []
end
println("Running solver and storing results...")

directSolverNonLinearBar(
	linear_problem,
	dataset;
)
for num_data_pts in 8:4:max_datapoints
	dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
	inner_results = []

	for num_ele in 2:1:max_elements
		@show num_data_pts, num_ele
		times = []
		costs = []
		if any(r -> r["Datapoints"] == num_data_pts && r["Elements"] == num_ele, results_list)
			continue
		end
		for repetitions in 1:100

			linear_problem, _ = get_problems(num_ele)
			t1 = time()
			result = directSolverNonLinearBar(
				linear_problem,
				dataset;
			)
			t2 = time()
			push!(times, t2 - t1)
			push!(costs, result.cost[end])
			# Log result
		end
		push!(inner_results, Dict(
			"Datapoints" => num_data_pts,
			"Elements" => num_ele,
			"Mean Time" => mean(times),
			"Mean cost" => mean(costs),
		))
	end
	# Save after each dataset sweep
	append!(results_list, inner_results)
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end


pivoted = unstack(DataFrame(results_list), :Elements, :Datapoints, "Mean Time")

# === Plot ===

elements = pivoted.Elements
datapoints = names(pivoted)[2:end]  # Skip :Elements
heatmap_data = Matrix(select(pivoted, Not(:Elements)))

heatmap(datapoints, elements, heatmap_data,
	xlabel = "Data points",
	ylabel = "Elements",
	colorbar_title = "Mean Time",
	c = :viridis,
)


savefig(replace(results_file, "results.json" => "heatmap.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))
