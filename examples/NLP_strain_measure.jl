include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables


results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	# Solve a problem to ensure precompilation
	NLP_solver(
		linear_problem,
		dataset;
		use_L1_norm = true,
		random_init_data = false,
		parameter_file = "NLP_params.prm")
	linear_problem, nonlinear_problem = get_problems(num_ele)
	for is_non_linear in [true, false]
		for random_init in [false, true]
			t1 = time()
			result = NLP_solver(
				is_non_linear ? nonlinear_problem : linear_problem,
				dataset;
				use_L1_norm = true,
				random_init_data = random_init,
				parameter_file = "NLP_params.prm",
			)
			t2 = time()

			push!(results_list, Dict(
				"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
				"Strain measure" => is_non_linear ? "Nonlinear" : "Linear",
				"Work" => result.solvetime[1], "Result" => result,
			))
		end
	end
	# Save results to file
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end
df = DataFrame(results_list)

unstack(df, "Strain measure", "Initialization", "Work")
# pivot(df, "Strain measure", "Initialization", "Work")

table = process_results(df, results_file, ("Work", "Work"))
table_colswapped = table[:, ["Strain measure", "No initialization", "Nullspace initialization"]]
table_rowswapped = table_colswapped[end:-1:1, :]


open(replace(results_file, ".json" => ".tex"), "w") do f
	pretty_table(f, table_rowswapped, backend = Val(:latex), show_subheader = false)
	pretty_table(table_rowswapped, show_subheader = false)
end


uncomment_pgfplotsset_blocks(dirname(results_file))
