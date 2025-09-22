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
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")

	linear_problem, nonlinear_problem = get_problems(num_ele)
	for is_L1 in [true, false]

		for random_init in [false, true]
			t1 = time()
			result = NLP_solver(
				linear_problem,
				dataset;
				use_L1_norm = is_L1,
				random_init_data = random_init,
				parameter_file = "NLP_params.prm")
			t2 = time()

			push!(results_list, Dict(
				"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
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
uncomment_pgfplotsset_blocks(dirname(results_file))


# # Convert results to DataFrame
# df = DataFrame(results_list)

# begin
# 	init_groups = groupby(df, ["Initialization"])
# 	p = plot(yscale = :log2, legend = :topright, xlabel = "Objective function", ylabel = "Solve time (s)")
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

# # res = SolveResults(; Dict(Symbol(k) => v for (k, v) in df[1, :Result])...)
# gr()
# plot_results(df[5, :Result], dataset = dataset)
# # process_results(df, results_file)

