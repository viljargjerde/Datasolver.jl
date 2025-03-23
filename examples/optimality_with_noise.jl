include("basic_setup.jl")
using Plots
using StatsPlots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables



joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
datasets = [create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit, noise_magnitude = 0.1) for _ in 1:100]
# Define shorthand functions for the solvers
NLP(lin_problem, random_init, data) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, data; use_L1_norm = false, random_init_data = random_init)
directSolver(lin_problem, random_init, data) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	# Solve a problem to ensure precompilation
	NLP(false, false, datasets[1])
	directSolver(false, false, datasets[1])
	for i in eachindex(datasets)
		for use_NLP_solver in [false, true]
			for lin_problem in [true, false]
				for random_init in [false, true]
					t1 = time()
					result = use_NLP_solver ? NLP(lin_problem, random_init, datasets[i]) : directSolver(lin_problem, random_init, datasets[i])
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

groups = groupby(df, ["Solver", "Initialization", "Problem"])
groups[1]

###* Linear ###
begin
	p = plot(yscale = :log2, legend = false)

	for group in groups
		if "Linear" in group[!, "Problem"]
			continue
		end
		name = "$(group[1,"Solver"]) $(group[1,"Initialization"])"
		@show name
		boxplot!([name], group[!, "Cost"], markersize = 3)
	end
	savefig(joinpath(@__DIR__, "boxplot.tex"))
	p
end

# ###* Non-linear ###
# p = plot(yscale = :log2, legend = false)
# for group in groups
# 	if "Linear" ∉ group[!, "Problem"]
# 		continue
# 	end
# 	name = "$(group[1,"Solver"]) $(group[1,"Initialization"])"
# 	@show name
# 	boxplot(repeat([name]), group[!, "Cost"])
# end
# p

# select!(df, ["Problem", "Initialization", "Solver", "Cost"])
# process_results_3vars(df, results_file, y = ("Cost", "Cost"), scatter = true)

# res = directSolver(true,false,datasets[1])
# plot_results(res,dataset=datasets[1])

