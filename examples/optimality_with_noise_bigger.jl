include("basic_setup.jl")
using Plots
using StatsPlots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables



results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []

datasets = [create_dataset(numDataPts, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = 0.01) for _ in 1:100]

linear_problem, nonlinear_problem = get_problems(4)
scatter(datasets[1].E, datasets[1].S)
# Define shorthand functions for the solvers
NLP(lin_problem, random_init, data) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, data; use_L1_norm = false, random_init_data = random_init, worklimit = 200.0, parameter_file = "NLP_params.prm")
directSolver(lin_problem, random_init, data) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init, NR_tol = 1e-9)
greedySearch(lin_problem, random_init, data) = Datasolver.greedyLocalSearchSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init, search_iters = 1000, NR_tol = 1e-9)

# solvers = ["ADM", "GO-ADM"]
solvers = ["ADM", "GO-ADM"]
num_exps = length(datasets) * length(solvers) * 2 * 2


if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	results_list = Vector{Dict}(undef, num_exps)
	Threads.@threads for i in eachindex(datasets)
		di = 0
		for lin_problem in [true, false]
			for random_init in [true, false]
				for solver in solvers
					@show solver
					t1 = time()
					if solver == "ADM"
						result = directSolver(lin_problem, random_init, datasets[i])
					elseif solver == "GO-ADM"
						result = greedySearch(lin_problem, random_init, datasets[i])
					elseif solver == "MINLP"
						result = NLP(lin_problem, random_init, datasets[i])
					end

					t2 = time()
					@show result.cost[end]
					results_list[length(solvers)*4*(i-1)+1+di] =
						Dict(
							"Solver" => solver,
							"Initialization" => random_init ? "Random" : "Nullspace",
							"Problem" => lin_problem ? "Linear" : "Nonlinear",
							"Cost" => result.cost[end],
							"Total solve time" => t2 - t1,
							"Model solve time" => result.solvetime[2],
							"Work" => result.solvetime[1],
							"Result" => result,
							"Dataset" => i,
							"S" => datasets[i].S,
							"E" => datasets[i].E,
						)
					di += 1
				end
			end
		end
	end
	# Save results to file
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end

df = DataFrame(Dict.(results_list))
combos = unique(select(df, [:Initialization, :Problem]))

function latex_sci(x::Real; digits = 1)
	if x == 0
		return "0"
	end
	exponent = floor(Int, log10(abs(x)))
	base = round(x / 10^exponent, digits = digits)
	return "\$$(base) \\times 10^{$(exponent)}\$"
end
# Create 4 subplots
plot_list = []
row_titles = []
col_titles = []

ymin = minimum(df.Cost)
ymax = maximum(df.Cost)
yticks = range(ymin, ymax; length = 5)


for row in eachrow(combos)
	init = row.Initialization
	prob = row.Problem
	push!(row_titles, "$(prob)")
	push!(col_titles, "$(init)")
	df_sub = filter(r -> r.Initialization == init && r.Problem == prob, df)
	p = @df df_sub boxplot(:Solver, :Cost, legend = false, yticks = yticks, ylim = (ymin, ymax), yformatter = latex_sci)
	println("$prob $init")
	println("ADM: $(sum(filter(r -> r.Solver == "ADM", df_sub).Cost))")
	println("GO-ADM: $(sum(filter(r -> r.Solver == "GO-ADM", df_sub).Cost))")
	println()
	push!(plot_list, p)
end


row_titles = unique(row_titles)
col_titles = unique(col_titles)

plt = plot(plot_list[1:4]..., layout = (2, 2), size = (400, 600))
title!(plt[1], col_titles[1], titlefont = font(12))
title!(plt[2], col_titles[2], titlefont = font(12))

# Add row titles (on y-axis)
ylabel!(plt[1], row_titles[1])
ylabel!(plt[3], row_titles[2])

savefig(plt, replace(results_file, "results.json" => "boxplot.tex"))



############# Table ###############
grouped = groupby(df, [:Initialization, :Problem, :Solver])

stats_df = combine(grouped,
	:Cost => mean => :mean,
	:Cost => median => :median,
	:Cost => std => :std,
	:Cost => minimum => :min,
	:Cost => maximum => :max,
)


table = select(filter(row -> row.Initialization == "Nullspace" && row.Problem == "Linear", stats_df), Not([:Initialization, :Problem]))
pretty_table(table, show_subheader = false)
open(replace(results_file, "results.json" => "nullspace-linear.tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
end

table = select(filter(row -> row.Initialization == "Nullspace" && row.Problem == "Nonlinear", stats_df), Not([:Initialization, :Problem]))
pretty_table(table, show_subheader = false)
open(replace(results_file, "results.json" => "nullspace-nonlinear.tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
end
uncomment_pgfplotsset_blocks(dirname(results_file))
