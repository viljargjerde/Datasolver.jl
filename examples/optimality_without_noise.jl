include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables



results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, nonlinear_problem = get_problems(4)

# Define shorthand functions for the solvers
NLP(lin_problem, random_init) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, dataset; use_L1_norm = false, random_init_data = random_init, parameter_file = "NLP_params.prm")
directSolver(lin_problem, random_init) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, dataset; random_init_data = random_init)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	for use_NLP_solver in [false, true]
		for lin_problem in [true, false]
			for random_init in [false, true]
				t1 = time()
				result = use_NLP_solver ? NLP(lin_problem, random_init) : directSolver(lin_problem, random_init)
				t2 = time()

				push!(
					results_list,
					Dict(
						"Solver" => use_NLP_solver ? "NLP" : "Direct",
						"Initialization" => random_init ? "Random" : "Nullspace",
						"Problem" => lin_problem ? "Linear" : "Nonlinear",
						"Cost" => result.cost[end],
						"Result" => result,
					),
				)
			end
		end
	end
	# Save results to file
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end

# Convert results to DataFrame
df = DataFrame(Dict.(results_list))
select!(df, ["Problem", "Initialization", "Solver", "Cost"])
df_direct = filter(row -> row.Solver == "Direct", df)
df_NLP = filter(row -> row.Solver == "NLP", df)

df_NLP = unstack(df_NLP, :Problem, :Initialization, :Cost)
df_direct = unstack(df_direct, :Problem, :Initialization, :Cost)

# gr()
scatter(df_direct[:, "Problem"], df_direct[:, "Nullspace"], marker = :star, label = "Direct solver, nullspace", scale = :log2, legend = :left)
scatter!(df_direct[:, "Problem"], df_direct[:, "Random"], marker = :square, label = "Direct solver, random")
scatter!(df_NLP[:, "Problem"], df_NLP[:, "Nullspace"], marker = :circle, label = "NLP solver, nullspace")
scatter!(df_NLP[:, "Problem"], df_NLP[:, "Random"], marker = :star, label = "NLP solver, random")

savefig(replace(results_file, "results.json" => "figure.tex"))
# process_results_3vars(df, results_file, y = ("Cost", "Cost"), scatter = true)

# plot(get_final(df_direct[1, :]["Result"]).u)
# plot!(get_final(df_direct[3, :]["Result"]).u)
# gr()
# plot_results(df_direct[3, :]["Result"], dataset = dataset)
# plot_dataset(dataset)

# df_direct
