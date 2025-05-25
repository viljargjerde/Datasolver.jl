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

directSolver(prob, data) = directSolverNonLinearBar(prob, data; NR_tol = 1e-8)
greedySearch(prob, data) = Datasolver.greedyLocalSearchSolverNonLinearBar(prob, data; search_iters = 100, NR_tol = 1e-8)

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
			open(results_file, "w") do f
				JSON.print(f, ratios)
			end
		end
	end
end


df = DataFrame(Dict.(ratios))
pivot = unstack(df, :Elements, :Datapoints, :Ratio)
x = parse.(Int, names(pivot)[2:end])  # Datapoints
y = pivot[:, 1]          # Elements
z = Matrix(pivot[:, 2:end])  # Ratios
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
