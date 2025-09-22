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
# create_dataset(16, x -> bar_E * x, -strain_limit, strain_limit)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")

	linear_problem, nonlinear_problem = get_problems(num_ele)
	for is_L1 in [false, true]

		for lin_problem in [true, false]
			t1 = time()
			result = NLP_solver(
				lin_problem ? linear_problem : nonlinear_problem,
				dataset;
				use_L1_norm = is_L1,
				random_init_data = false,
			)
			t2 = time()

			push!(results_list, Dict(
				"Strain measure" => lin_problem ? "Linear" : "Non linear",
				"Objective function" => is_L1 ? "L1" : "L2",
				# "Solve time" => t2 - t1,
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
# process_results(df, results_file, ("Work", "Work"))

# begin
# 	init_groups = groupby(df, ["Initialization"])
# 	p = plot(yscale = :log2, legend = :topright, xlabel = "Objective function", ylabel = "Work")
# 	for init_group in init_groups
# 		initialization = init_group[1, "Initialization"]
# 		L1 = median(init_group[init_group[:, "Objective function"].=="L1", "Solve time"])
# 		L2 = median(init_group[init_group[:, "Objective function"].=="L2", "Solve time"])
# 		plot!(["L1", "L2"], [L1, L2], label = initialization, marker = :circle, markercolor = :transparent)
# 		@show L1
# 	end
# 	savefig(replace(results_file, "results.json" => "lineplot.tex"))
# 	p
# end

# res = SolveResults(; Dict(Symbol(k) => v for (k, v) in df[1, :Result])...)
# gr()
# plot_results(df[5, :Result], dataset = dataset)
# process_results(df, results_file)

