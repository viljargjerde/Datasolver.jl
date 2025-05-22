include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using Interpolations, StatsBase
using LaTeXStrings
using ProgressBars

num_data_pts = 2^10
num_ele = 2^8

dataset = create_dataset(num_data_pts, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = 0.01)

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, nonlinear_problem = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	# Solve a problem to ensure precompilation
	Datasolver.greedyLocalSearchSolverNonLinearBar(
		nonlinear_problem,
		dataset,
		random_init_data = false,
		NR_tol = 1e-8,
	)
	@show num_ele
	for i in ProgressBar(1:100)
		dataset = create_dataset(num_data_pts, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = 0.01)
		t1 = time()
		result = greedyLocalSearchSolverNonLinearBar(
			nonlinear_problem,
			dataset;
			random_init_data = false,
			NR_tol = 1e-8,
			cache_ADM = true,
			search_iters = 200,
		)
		t2 = time()

		push!(results_list, Dict(
			"Elements" => num_ele,
			"Solve time" => t2 - t1,
			"Result" => result,
		))
	end
	# Save results to file
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end

# Convert results to DataFrame
df = DataFrame(results_list)

# process_results(df, results_file)


begin

	p = plot(xlabel = "Time (s)", ylabel = "Cost")
	# plot(results_list[1]["Result"].solvetime, results_list[1]["Result"].cost, label = nothing)
	for res in results_list
		plot!(res["Result"]["solvetime"][begin:end-1], res["Result"]["cost"], label = nothing, alpha = 0.3)
	end

	##### Mean #######
	t_max = maximum([res["Result"]["solvetime"][end-1] for res in results_list])
	t_final = mean([res["Result"]["solvetime"][end] for res in results_list])
	# Define a common time grid
	t_common = range(0, stop = t_max, length = 200)

	# Interpolate each result onto the common time grid
	cost_matrix = []

	for res in results_list
		t = res["Result"]["solvetime"][begin:end-1]
		c = res["Result"]["cost"]

		if length(t) > 1  # skip if data is too sparse
			itp = LinearInterpolation(t, c, extrapolation_bc = Flat())
			push!(cost_matrix, itp.(t_common))
		end
	end

	# Convert to matrix
	cost_array = hcat(cost_matrix...)'

	# Compute mean or median across runs
	mean_cost = mean(cost_array, dims = 1)
	median_cost = mapslices(median, cost_array; dims = 1)

	# Plot the mean or median on top
	plot!(t_common, vec(mean_cost), label = "Mean", lw = 2, color = :black, legend = :topright)
	# Or for median:
	# plot!(t_common, vec(median_cost), label = "Median", lw = 2, color = :blue)
	mask = t_common .> 0.05
	t_fit = t_common[mask]
	y_fit = vec(mean_cost[mask])
	a, b, f = estimate_powerlaw(t_fit, y_fit)
	plot!(t_fit, f.(t_fit), lw = 1, label = L"y = %$(round(Int,a))t^{%$(round(b, digits = 2))}", color = :red, linestyle = :dash)
	savefig(replace(results_file, "results.json" => "figure.tex"))
	update_tex_command(all_results_file, "SAADMpowerlaw", format(FormatExpr("y = {}t^{{{:.2f}}}"), round(Int, a), b))
	update_tex_command(all_results_file, "SAADMMeanFinalTime", format(FormatExpr("y = {:.1f}"), t_final))

	p
end

uncomment_pgfplotsset_blocks(dirname(results_file))

