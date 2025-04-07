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
# datasets = [create_dataset(numDataPts, x -> bar_E * tanh.(20x), -strain_limit, strain_limit, noise_magnitude = 0.2) for _ in 1:10]

linear_problem, nonlinear_problem = get_problems(4)
scatter(datasets[1].E, datasets[1].S)
# Define shorthand functions for the solvers
NLP(lin_problem, random_init, data) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, data; use_L1_norm = false, random_init_data = random_init, worklimit = 200.0, parameter_file = "NLP_params.prm")
directSolver(lin_problem, random_init, data) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init)
localSearch(lin_problem, random_init, data) = Datasolver.localSearchSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init, search_iters = 1000, neighborhood_size = 2)
greedySearch(lin_problem, random_init, data) = Datasolver.greedyLocalSearchSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init, max_search_iters = 1000)
hybridSearch(lin_problem, random_init, data) = Datasolver.hybridLocalSearchSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init, search_iters = 1000, greedy_prob = 0.9)

solvers = ["Direct", "Greedy"]
# solvers = ["Direct", "Greedy", "NLP"]
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
					if solver == "Direct"
						result = directSolver(lin_problem, random_init, datasets[i])
					elseif solver == "Search"
						result = localSearch(lin_problem, random_init, datasets[i])
					elseif solver == "Greedy"
						result = greedySearch(lin_problem, random_init, datasets[i])
					elseif solver == "Hybrid"
						result = hybridSearch(lin_problem, random_init, datasets[i])
					elseif solver == "NLP"
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
	println("Direct: $(sum(filter(r -> r.Solver == "Direct", df_sub).Cost))")
	println("Search: $(sum(filter(r -> r.Solver == "Search", df_sub).Cost))")
	println("Greedy: $(sum(filter(r -> r.Solver == "Greedy", df_sub).Cost))")
	println("NLP: $(sum(filter(r -> r.Solver == "NLP", df_sub).Cost))")
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

#############       ###############

df
df_direct = filter(row -> row.Solver == "Direct" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
df_search = filter(row -> row.Solver == "Search" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
df_greedy = filter(row -> row.Solver == "Greedy" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
df_NLP = filter(row -> row.Solver == "NLP" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
sort!(df_direct, :Dataset)
sort!(df_search, :Dataset)
sort!(df_greedy, :Dataset)
sort!(df_NLP, :Dataset)
cost_diff = abs.(df_direct.Cost .- df_greedy.Cost)
# cost_diff = abs.(df_direct.Cost .- df_NLP.Cost)
max_diff_idx = argmax(cost_diff)
DE = Vector{Float64}(df_direct[max_diff_idx, :E])
DS = Vector{Float64}(df_direct[max_diff_idx, :S])
scatter(DE, DS, label = "Data", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, line = :dash, markercolor = paired_colors[1])

s = Datasolver.get_initialization_s(linear_problem)
best_idxs = Datasolver.find_closest_idx(DS, s)
S = DS[best_idxs]
E = DE[best_idxs]
scatter!(E, S, label = "Initialization", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, markercolor = paired_colors[5], line = :dash)

r1 = df_direct[max_diff_idx, :Result]
scatter!(r1["e"][end], r1["s"][end], label = "Direct", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :diamond, markercolor = paired_colors[9], line = :dash)

r3 = df_greedy[max_diff_idx, :Result]
scatter!(r3["e"][end], r3["s"][end], label = "Greedy", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :circle, markercolor = paired_colors[12], line = :dash)

r2 = df_NLP[max_diff_idx, :Result]
scatter!(r2["e"][end], r2["s"][end], label = "NLP", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :star, markersize = 2, markercolor = paired_colors[9], line = :dash)

savefig(replace(results_file, "results.json" => "example-solution.tex"))

# r["e"][end] .== r["E"][end]

# plot_results(SolveResults(; Dict(Symbol(k) => v for (k, v) in df_direct[max_diff_idx, :Result])...), dataset = datasets[max_diff_idx])


# Combine into 2x2 grid

# plot_dataset(datasets[2])
# plot(plot_list[1])
# plot(plot_list[2])
# plot(plot_list[3])
# plot(plot_list[4])


# df.label = string.(df.Solver, " | ", df.Initialization, " | ", df.Problem)

# @df df boxplot(:label, :Cost, xlabel = "Solver | Initialization | Problem", ylabel = "Cost", legend = :topright, scale = :log2)

# gr()
# r_NLP = filter(row -> row.Solver == "NLP", df)
# r_direct = filter(row -> row.Solver == "Direct", df)
# r_search = filter(row -> row.Solver == "Search", df)
# @df r_NLP boxplot(:Problem, :Cost, group = :Initialization, label = "NLP", color = :blue)
# @df r_direct boxplot!(:Problem, :Cost, group = :Initialization, label = "Direct", color = :orange)
# @df r_search boxplot!(:Problem, :Cost, group = :Initialization, label = "Search", color = :green)


# r = combine(groupby(df, ["Initialization", "Solver", "Problem"]), :Cost => mean)
# r_NLP = select(filter(row -> row.Solver == "NLP", r), Not(:Solver))
# r_direct = select(filter(row -> row.Solver == "Direct", r), Not(:Solver))
# sort!(r_NLP, :Problem)
# sort!(r_direct, :Problem)
# r_NLP = unstack(r_NLP, :Problem, :Initialization, :Cost_mean)
# r_direct = unstack(r_direct, :Problem, :Initialization, :Cost_mean)

# select!(df, ["Problem", "Initialization", "Solver", "Cost"])
# df_direct = filter(row -> row.Solver == "Direct", df)
# df_NLP = filter(row -> row.Solver == "NLP", df)

# df_NLP = unstack(df_NLP, :Problem, :Initialization, :Cost)
# df_direct = unstack(df_direct, :Problem, :Initialization, :Cost)

# scatter(df_direct[:, "Problem"], df_direct[:, "Nullspace"], marker = :star, label = "Direct solver, nullspace", scale = :log2, legend = :left)
# scatter!(df_direct[:, "Problem"], df_direct[:, "Random"], marker = :square, label = "Direct solver, random")
# scatter!(df_NLP[:, "Problem"], df_NLP[:, "Nullspace"], marker = :circle, label = "NLP solver, nullspace")
# scatter!(df_NLP[:, "Problem"], df_NLP[:, "Random"], marker = :star, label = "NLP solver, random")

# savefig(replace(results_file, ".json" => ".tex"))

# # Convert results to DataFrame
# df = DataFrame(Dict.(results_list))

# groups = groupby(df, ["Solver", "Initialization", "Problem"])

# ###* Linear ###
# begin
# 	p = plot(yscale = :log2, legend = false)

# 	for group in groups
# 		if "Linear" in group[!, "Problem"]
# 			continue
# 		end
# 		name = "$(group[1,"Solver"]) $(group[1,"Initialization"])"
# 		@show name
# 		boxplot!([name], group[!, "Cost"], markersize = 3)
# 	end
# 	savefig(joinpath(@__DIR__, "boxplot.tex"))
# 	p
# end

# ###* Nonlinear ###
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

