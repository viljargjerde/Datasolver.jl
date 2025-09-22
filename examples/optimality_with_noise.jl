include("basic_setup.jl")
using Plots
using StatsPlots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using Printf



results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
datasets = [create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit, noise_magnitude = 0.1) for _ in 1:100]
linear_problem, nonlinear_problem = get_problems(4)

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

					push!(
						results_list,
						Dict(
							"Solver" => use_NLP_solver ? "NLP" : "Direct",
							"Initialization" => random_init ? "Random" : "Nullspace",
							"Problem" => lin_problem ? "Linear" : "Non-linear",
							"Cost" => result.cost[end],
							"Total solve time" => t2 - t1,
							"Model solve time" => result.solvetime[2],
							"Work" => result.solvetime[1],
							"Result" => result,
							"Dataset" => i,
							"S" => datasets[i].S,
							"E" => datasets[i].E,
						),
					)
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

# Create 4 subplots
plot_list = []
row_titles = []
col_titles = []

ymin = minimum(df.Cost)
ymax = maximum(df.Cost)
yticks = range(ymin, ymax; length = 5)

function latex_sci(x::Real; digits = 1)
	if x == 0
		return "0"
	end
	exponent = floor(Int, log10(abs(x)))
	base = round(x / 10^exponent, digits = digits)
	return "\$$(base) \\times 10^{$(exponent)}\$"
end

for row in eachrow(combos)
	init = row.Initialization
	prob = row.Problem
	push!(row_titles, "$(prob)")
	push!(col_titles, "$(init)")
	df_sub = filter(r -> r.Initialization == init && r.Problem == prob, df)

	p = @df df_sub boxplot(:Solver, :Cost, legend = false, yticks = yticks, ylim = (ymin, ymax), yformatter = latex_sci)


	push!(plot_list, p)
end
row_titles = unique(row_titles)
col_titles = unique(col_titles)

plt = plot(plot_list[1:4]..., layout = (2, 2), size = (400, 600))
title!(plt[1], col_titles[1], titlefontsize = 10, titlefont = font(12))
title!(plt[2], col_titles[2], titlefont = font(12))

# Add row titles (on y-axis)
ylabel!(plt[1], row_titles[1])
ylabel!(plt[3], row_titles[2])

savefig(plt, replace(results_file, "results.json" => "boxplot.tex"))


df
df_direct = filter(row -> row.Solver == "Direct" && row.Problem == "Non-linear" && row.Initialization == "Nullspace", df)
df_NLP = filter(row -> row.Solver == "NLP" && row.Problem == "Non-linear" && row.Initialization == "Nullspace", df)
sort!(df_direct, :Dataset)
sort!(df_NLP, :Dataset)
cost_diff = abs.(df_direct.Cost .- df_NLP.Cost)
max_diff_idx = argmax(cost_diff)
df_direct[max_diff_idx, :]
df_NLP[max_diff_idx, :]
DE = Vector{Float64}(df_direct[max_diff_idx, :E])
DS = Vector{Float64}(df_direct[max_diff_idx, :S])
gr()
scatter(DE, DS, label = "Data", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, line = :dash)

s = Datasolver.get_initialization_s(linear_problem)
best_idxs = Datasolver.find_closest_idx(DS, s)
S = DS[best_idxs]
E = DE[best_idxs]
scatter!(E, S, label = "Initialization", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, markercolor = :yellow, line = :dash)

r1 = df_direct[max_diff_idx, :Result]
scatter!(r1["e"][end], r1["s"][end], label = "Direct", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :circle, markercolor = :red, line = :dash)
@show sum(r1["s"][end])
@show sum(r1["S"][end])
@show r1["S"][end] .- (r1["s"][end])
@show sum(r1["S"][end] .- (r1["s"][end]))

r2 = df_NLP[max_diff_idx, :Result]
@show sum(r2["s"][end])
@show sum(r2["S"][end])
@show sum(r2["S"][end] .- (r2["s"][end]))
@show r2["S"][end] .- (r2["s"][end])

@show r1["s"][end] .- r2["s"][end]
@show r1["S"][end] .- r2["S"][end]

scatter!(r2["e"][end], r2["s"][end], label = "NLP", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :star, markercolor = :green, line = :dash)


# r["e"][end] .== r["E"][end]

# plot_results(SolveResults(; Dict(Symbol(k) => v for (k, v) in df_direct[max_diff_idx, :Result])...), dataset = datasets[max_diff_idx])


# Combine into 2x2 grid

# # plot_dataset(datasets[2])
# plot(plot_list[1])
# plot(plot_list[2])
# plot(plot_list[3])
# plot(plot_list[4])


# df.label = string.(df.Solver, " | ", df.Initialization, " | ", df.Problem)

# @df df boxplot(:label, :Cost, xlabel = "Solver | Initialization | Problem", ylabel = "Cost", legend = :topright,scale=:log2)

# gr()
# r_NLP = filter(row -> row.Solver == "NLP", df)
# r_direct = filter(row -> row.Solver == "Direct", df)
# @df r_NLP boxplot(:Problem, :Cost, group = :Initialization, label = "NLP", color = :blue)
# @df r_direct boxplot!(:Problem, :Cost, group = :Initialization, label = "Direct", color = :orange)


# r  = combine(groupby(df,["Initialization","Solver","Problem"]), :Cost => mean)
# r_NLP = select(filter(row -> row.Solver == "NLP", r),Not(:Solver))
# r_direct = select(filter(row -> row.Solver == "Direct", r),Not(:Solver))
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

