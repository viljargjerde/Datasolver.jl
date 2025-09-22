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

# datasets = [create_dataset(128, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = 0.01) for _ in 1:100]
# datasets = [create_dataset(numDataPts, x -> bar_E * tanh.(20x), -strain_limit, strain_limit, noise_magnitude = 0.2) for _ in 1:10]

directSolver(prob, data) = directSolverNonLinearBar(prob, data; NR_tol = 1e-8)
greedySearch(prob, data) = Datasolver.greedyLocalSearchSolverNonLinearBar(prob, data; search_iters = 100, NR_tol = 1e-8)

# solvers = ["ADM", "GO-ADM"]
solvers = ["ADM", "GO-ADM"]

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")

num_ele_list = [2^n for n in 2:10]
numDataPts_list = [2^n for n in 2:10]
num_reps = 100

# ratios = Dict()
ratios = []
if isfile(results_file)
	println("Loading existing results from file...")
	ratios = JSON.parsefile(results_file)
else
	for num_ele in num_ele_list
		for numDataPts in numDataPts_list
			println("Running for num_ele = $num_ele, numDataPts = $numDataPts")
			linear_problem, nonlinear_problem = get_problems(num_ele)

			local_ratios = zeros(num_reps)

			Threads.@threads for rep in 1:num_reps
				data = create_dataset(numDataPts, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = 0.01)
				try

					# Direct solver
					direct_result = directSolver(nonlinear_problem, data)

					# Greedy solver
					greedy_result = greedySearch(nonlinear_problem, data)

					local_ratios[rep] = direct_result.cost[end] / greedy_result.cost[end]
				catch e
					local_ratios[rep] = 1
					println("Failed for one experiment with num_ele = $num_ele, numDataPts = $numDataPts, skipping $rep")
				end
			end

			avg_ratio = mean(local_ratios)
			push!(ratios, Dict(
				"Elements" => num_ele,
				"Datapoints" => numDataPts,
				"Ratio" => avg_ratio,
			))
			# ratios[(num_ele, numDataPts)] = avg_ratio
			open(results_file, "w") do f
				JSON.print(f, ratios)
			end
		end
	end
end
# Save the ratios for plotting later


df = DataFrame(Dict.(ratios))
pivot = unstack(df, :Elements, :Datapoints, :Ratio)
x = parse.(Int, names(pivot)[2:end])  # Datapoints
y = pivot[:, 1]          # Elements
z = Matrix(pivot[:, 2:end])  # Ratios
# z[1, 3] = NaN
heatmap(
	x[2:end],
	y[2:end],
	z[2:end, :1:end-1], ;
	xticks = [x[i] for i in 1:2:length(x)],
	yticks = [y[i] for i in 1:2:length(y)],
	xlabel = "Data points",
	ylabel = "Elements",
	colorbar_title = "ADM / Greedy Cost Ratio",
	scale = :log2,
	c = :viridis,
)

savefig(replace(results_file, "results.json" => "heatmap.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))

# combos = unique(select(df, [:Initialization, :Problem]))

# function latex_sci(x::Real; digits = 1)
# 	if x == 0
# 		return "0"
# 	end
# 	exponent = floor(Int, log10(abs(x)))
# 	base = round(x / 10^exponent, digits = digits)
# 	return "\$$(base) \\times 10^{$(exponent)}\$"
# end
# # Create 4 subplots
# plot_list = []
# row_titles = []
# col_titles = []

# ymin = minimum(df.Cost)
# ymax = maximum(df.Cost)
# yticks = range(ymin, ymax; length = 5)


# for row in eachrow(combos)
# 	init = row.Initialization
# 	prob = row.Problem
# 	push!(row_titles, "$(prob)")
# 	push!(col_titles, "$(init)")
# 	df_sub = filter(r -> r.Initialization == init && r.Problem == prob, df)
# 	p = @df df_sub boxplot(:Solver, :Cost, legend = false, yticks = yticks, ylim = (ymin, ymax), yformatter = latex_sci)
# 	println("$prob $init")
# 	println()
# 	push!(plot_list, p)
# end


# row_titles = unique(row_titles)
# col_titles = unique(col_titles)

# plt = plot(plot_list[1:4]..., layout = (2, 2), size = (400, 600))
# title!(plt[1], col_titles[1], titlefont = font(12))
# title!(plt[2], col_titles[2], titlefont = font(12))

# # Add row titles (on y-axis)
# ylabel!(plt[1], row_titles[1])
# ylabel!(plt[3], row_titles[2])

# savefig(plt, replace(results_file, "results.json" => "boxplot.tex"))



# ############# Table ###############
# grouped = groupby(df, [:Initialization, :Problem, :Solver])

# stats_df = combine(grouped,
# 	:Cost => mean => :mean,
# 	:Cost => median => :median,
# 	:Cost => std => :std,
# 	:Cost => minimum => :min,
# 	:Cost => maximum => :max,
# )


# table = select(filter(row -> row.Initialization == "Nullspace" && row.Problem == "Linear", stats_df), Not([:Initialization, :Problem]))
# pretty_table(table, show_subheader = false)
# open(replace(results_file, "results.json" => "nullspace-linear.tex"), "w") do f
# 	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
# end

# table = select(filter(row -> row.Initialization == "Nullspace" && row.Problem == "Nonlinear", stats_df), Not([:Initialization, :Problem]))
# pretty_table(table, show_subheader = false)
# open(replace(results_file, "results.json" => "nullspace-nonlinear.tex"), "w") do f
# 	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
# end

#############       ###############

# df
# df_direct = filter(row -> row.Solver == "ADM" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
# df_greedy = filter(row -> row.Solver == "GO-ADM" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
# sort!(df_direct, :Dataset)
# sort!(df_greedy, :Dataset)
# cost_diff = abs.(df_direct.Cost .- df_greedy.Cost)
# # cost_diff = abs.(df_direct.Cost .- df_NLP.Cost)
# max_diff_idx = argmax(cost_diff)
# DE = Vector{Float64}(df_direct[max_diff_idx, :E])
# DS = Vector{Float64}(df_direct[max_diff_idx, :S])

# ### For dataset plot ### 
# scatter(DE, DS, label = nothing, xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, line = :dash, markercolor = paired_colors[1])
# savefig(replace(results_file, "results.json" => "dataset.tex"))
# #########################
# scatter(DE, DS, label = "Data", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, line = :dash, markercolor = paired_colors[1])
# s = Datasolver.get_initialization_s(linear_problem)
# best_idxs = Datasolver.find_closest_idx(DS, s)
# S = DS[best_idxs]
# E = DE[best_idxs]
# scatter!(E, S, label = "Initialization", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, markercolor = paired_colors[5], line = :dash)

# r1 = df_direct[max_diff_idx, :Result]
# scatter!(r1["e"][end], r1["s"][end], label = "ADM", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :diamond, markercolor = paired_colors[9], line = :dash)

# r2 = df_greedy[max_diff_idx, :Result]
# scatter!(r2["e"][end], r2["s"][end], label = "GO-ADM", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :circle, markercolor = paired_colors[12], line = :dash)

# savefig(replace(results_file, "results.json" => "example-solution.tex"))

