include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using JSON
using PrettyTables

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")

linear_problem, nonlinear_problem = get_problems(4)

NLP(lin_problem, random_init) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, dataset; use_L1_norm = false, random_init_data = random_init, parameter_file = "NLP_params.prm")
directSolver(lin_problem, random_init) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, dataset; random_init_data = random_init)
greedySolver(lin_problem, random_init) = Datasolver.greedyLocalSearchSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, dataset; random_init_data = random_init, search_iters = 1000)

solvers = ["MINLP", "ADM", "GO-ADM"]
problems = [true, false]  # true = Linear, false = Nonlinear
inits = [false, true]     # false = Nullspace, true = Random

n_random_repeats = 100

results_list = []

if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solvers and storing results...")
	for solver in solvers, lin_problem in problems, random_init in inits
		n_repeats = random_init ? n_random_repeats : 1

		for _ in 1:n_repeats
			t1 = time()
			result = solver == "MINLP" ? NLP(lin_problem, random_init) :
					 solver == "ADM" ? directSolver(lin_problem, random_init) :
					 greedySolver(lin_problem, random_init)
			t2 = time()

			push!(results_list, Dict(
				"Solver" => solver,
				"Initialization" => random_init ? "Random" : "Nullspace",
				"Problem" => lin_problem ? "Linear" : "Nonlinear",
				"Cost" => result.cost[end],
			))
		end
	end

	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end

using StatsPlots
df = DataFrame(results_list)
df_summary = combine(groupby(df, ["Solver", "Initialization", "Problem"]),
	:Cost => mean => :MeanCost)


df_linear = filter(:Problem => ==("Linear"), df_summary)

df_nonlinear = filter(:Problem => ==("Nonlinear"), df_summary)

ymin = minimum(df_summary.MeanCost)
ymax = maximum(df_summary.MeanCost)

p1 = @df df_linear groupedbar(
	:Solver, :MeanCost,
	group = :Initialization,
	bar_position = :dodge,
	legend = :top,
	ylabel = "Mean Cost",
	title = "Problem: Linear",
	ylim = (ymin, ymax),
)

p2 = @df df_nonlinear groupedbar(
	:Solver, :MeanCost,
	group = :Initialization,
	bar_position = :dodge,
	legend = nothing,
	ylabel = "Mean Cost",
	title = "Problem: Nonlinear",
	ylim = (ymin, ymax),
)

plot(p1, p2, layout = (2, 1), size = (400, 600))

savefig(replace(results_file, "results.json" => "figure.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))
