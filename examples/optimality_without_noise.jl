include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables



joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []

# Define shorthand functions for the solvers
NLP(lin_problem, random_init) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, dataset; use_L1_norm = false, random_init_data = random_init)
directSolver(lin_problem, random_init) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, dataset; random_init_data = random_init)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	# Solve a problem to ensure precompilation
	NLP(false, false)
	directSolver(false, false)
	for i in 1:10
		for use_NLP_solver in [false, true]
			for lin_problem in [true, false]
				for random_init in [false, true]
					t1 = time()
					result = use_NLP_solver ? NLP(lin_problem, random_init) : directSolver(lin_problem, random_init)
					t2 = time()

					push!(results_list, Dict(
						"Solver" => use_NLP_solver ? "NLP" : "Direct",
						"Initialization" => random_init ? "Random" : "Nullspace",
						"Problem" => lin_problem ? "Linear" : "Non-linear",
						"Cost" => result.cost[end],
						"Result" => result,
					))
				end
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

process_results_3vars(df, results_file, y = ("Cost", "Cost"), scatter = true)


