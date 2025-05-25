include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using StatsPlots
results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
num_ele = 5
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")

	linear_problem, nonlinear_problem = get_problems(num_ele)
	for is_L1 in [true, false]

		for lin_problem in [true, false]
			t1 = time()
			result = NLP_solver(
				lin_problem ? linear_problem : nonlinear_problem,
				dataset;
				use_L1_norm = is_L1,
				random_init_data = true,
				parameter_file = "NLP_params.prm",
				timelimit = 3600 * 5)
			t2 = time()

			push!(results_list, Dict(
				"Strain measure" => lin_problem ? "Linear" : "Nonlinear",
				"Objective function" => is_L1 ? "L1" : "L2",
				"Work" => result.solvetime[1],
				"Result" => result,
			))
		end
	end
	# Save results to file
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end


# # Convert results to DataFrame
begin
	df = DataFrame(results_list)

	table = unstack(select(df, Not(:Result)), "Strain measure", "Objective function", :Work)
	select!(table, ["Strain measure", "L1", "L2"])
	mat = Matrix(table[:, 2:end])

	# Labels
	row_labels = table[!, "Strain measure"]
	col_labels = names(table)[2:end]

	# Plot
	p = heatmap(col_labels, row_labels, log2.(mat), xlabel = "Objective function", ylabel = "Strain Measure", c = :viridis, colorbar_title = "Work (log2)")
	savefig(replace(results_file, "results.json" => "heatmap.tex"))
	p
end


open(replace(results_file, ".json" => ".tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
	pretty_table(table, show_subheader = false)
end

uncomment_pgfplotsset_blocks(dirname(results_file))
