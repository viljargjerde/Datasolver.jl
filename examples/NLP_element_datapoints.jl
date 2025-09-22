include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using StatsBase: Histogram, fit



results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")

max_elements = 33
max_datapoints = 256

if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	results_list = []
end
println("Running solver and storing results...")

for num_data_pts in 8:4:256
	dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
	timeout_occurred = false
	inner_results = []
	concecutive_timeouts = 0

	for num_ele in 2:1:max_elements

		# Check if this combination has already been solved
		if any(r -> r["Datapoints"] == num_data_pts && r["Elements"] == num_ele, results_list)
			continue
		end

		linear_problem, _ = get_problems(num_ele)
		t1 = time()
		result = NLP_solver(
			linear_problem,
			dataset;
			use_L1_norm = true,
			random_init_data = false,
			parameter_file = "NLP_params.prm",
			worklimit = 200,
		)
		t2 = time()

		# Check for timeout flag
		timeout_occurred = (result.solvetime[3] == 0.0)

		# Log result
		push!(inner_results, Dict(
			"Datapoints" => num_data_pts,
			"Elements" => num_ele,
			"Work" => result.solvetime[1],
			"Result" => result,
			"Timeout" => timeout_occurred,
		))
		if timeout_occurred
			concecutive_timeouts += 1
		else
			concecutive_timeouts = 0
		end

		# If timeout occurred, fill the rest of this dataset's loop with this result
		if concecutive_timeouts >= 3
			for fill_ele in (num_ele+1):1:max_elements
				result = deepcopy(result)
				push!(inner_results, Dict(
					"Datapoints" => num_data_pts,
					"Elements" => fill_ele,
					"Work" => result.solvetime[1],
					"Result" => result,
					"Timeout" => timeout_occurred,
				))
			end
			break
		end
	end

	# Save after each dataset sweep
	append!(results_list, inner_results)
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end


# Extract grid
data_pts_set = sort(unique([r["Datapoints"] for r in results_list]))
num_ele_set = sort(unique([r["Elements"] for r in results_list]))

# Prepare data matrices
Z = fill(NaN, length(data_pts_set), length(num_ele_set))           # Solve time
timeout_mask = falses(size(Z))                                     # Direct timeouts

# Fill matrices
for r in results_list
	i = findfirst(==(r["Datapoints"]), data_pts_set)
	j = findfirst(==(r["Elements"]), num_ele_set)
	Z[i, j] = r["Work"]
	timeout_mask[i, j] = r["Timeout"]
end


# === Plot ===


heatmap(
	num_ele_set, data_pts_set, Z;
	xlabel = "Number of Elements",
	ylabel = "Number of Data Points",
	colorbar_title = "Work Units",
	c = :viridis,
)

savefig(replace(results_file, "results.json" => "heatmap.tex"))
# timeout_mask = float.(timeout_mask)
# timeout_mask[timeout_mask.==false] .= NaN

# # --- Cliff overlay ---
# cliff_cmap = cgrad([RGBA(1, 0, 0, 1.0), RGBA(1, 0, 0, 1.0)])

# heatmap!(
# 	num_ele_set, data_pts_set, timeout_mask;
# 	color = cliff_cmap,
# 	colorbar = false,
# )
for (i, row) in enumerate(eachrow(timeout_mask))
	for j in eachindex(row)
		if row[j]
			println("Timeout at (", data_pts_set[i], ", ", num_ele_set[j], ") ratio: ", data_pts_set[i] / num_ele_set[j])
		end
	end
end


ratio = [r["Datapoints"] / r["Elements"] for r in results_list]
works = [r["Work"] for r in results_list]
scatter(ratio, works, markersize = 2, xlabel = "Ratio of datapoints to elements", ylabel = "Work", legend = false)
savefig(replace(results_file, "results.json" => "ratio_scatter.tex"))

ratio = ratio[works.>=200]
works = works[works.>=200]
histogram(ratio, bins = 20, xlabel = "Ratio of datapoints to elements", ylabel = "Number timed out", legend = false)
savefig(replace(results_file, "results.json" => "ratio_histogram.tex"))
###############



hist = StatsBase.fit(StatsBase.Histogram, ratio, nbins = 20)

for (e, w) in zip(hist.edges[1], hist.weights)
	if w > 0.0
		println("Ratio: ", e, " Count: ", w)
	end
end

# total_work = 0.0
# total_time = 0.0
# for res in df.Result
# 	work = res["solvetime"][1]
# 	solvetime = res["solvetime"][2]
# 	if work > 0.1 && solvetime > 0.1
# 		global total_work += work
# 		global total_time += solvetime
# 	end
# end

# println("Work to solvetime ratio: ", total_work / total_time)

# table = unstack(select(df, Not(["Result"])), :Initialization, :Work)


# xs_fitted = table[3, "Elements"]:table[end, "Elements"]
# p = plot(scale = :log2, legend = :topleft, xlabel = "Number of elements", ylabel = "Work units", palette = paired_colors) # :Paired_12 ,:tableau_20
# x_ticks = num_eles[3:2:end]
# plot!(table[3:end, "Elements"], table[3:end, "Nullspace initialization"], marker = :circle, label = "Nullspace initialization", xticks = x_ticks)
# a_null, b_null, f1 = estimate_powerlaw(table[3:end-1, "Elements"], table[3:end-1, "Nullspace initialization"])
# # a_null, b_null, c_null, f1 = estimate_quadratic_powerlaw(table[3:end, "Elements"], table[3:end, "Nullspace initialization"])

# update_tex_command(all_results_file, "NLPElementsPowerlawNull", format(FormatExpr("y \\propto x^{{{:.2f}}}"), b_null))
# update_tex_command(all_results_file, "NLPElementsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))
# # update_tex_command(all_results_file, "NLPElementsPowerlawNull", format(FormatExpr("y = {:.2g}x^{{{:.2f}}}"), a_null, b_null))

# plot!(xs_fitted, f1.(xs_fitted), label = nothing, linestyle = :dash)
# a_rand, b_rand, f2 = estimate_powerlaw(table[3:end-1, "Elements"], table[3:end-1, "Random initialization"])
# # a_rand, b_rand, c_rand, f2 = estimate_quadratic_powerlaw(table[3:end, "Elements"], table[3:end, "Random initialization"])

# update_tex_command(all_results_file, "NLPElementsPowerlawRand", format(FormatExpr("y \\propto x^{{{:.2f}}}"), b_rand))
# update_tex_command(all_results_file, "NLPElementsSpeedRatio", format(FormatExpr("{:.1f}"), a_rand / a_null))
# # update_tex_command(all_results_file, "NLPDatapointsPowerlawRand", format(FormatExpr("y = {:.2f}x^{{{:.2f}}}"), a_rand, b_rand))

# plot!(table[3:end, "Elements"], table[3:end, "Random initialization"], marker = :circle, label = "Random initialization")
# plot!(xs_fitted, f2.(xs_fitted), label = nothing, linestyle = :dash)
# # plot!(table[3:end, "Elements"], a_rand .* table[3:end, "Elements"] .^ b_rand, label = nothing, linestyle = :dash)
# savefig(replace(results_file, "results.json" => "lineplot.tex"))
# p









# # Convert results to DataFrame
# df = DataFrame(results_list)

# table = process_results(df, results_file)
# a_null, b_null, f1 = estimate_powerlaw(table[3:end, "Elements"], table[3:end, "Nullspace initialization"])
# a_rand, b_rand, f2 = estimate_powerlaw(table[3:end, "Elements"], table[3:end, "Random initialization"])
# f1(12)
# f1(200)
# f2(12)
# f2(200)
# b_rand
# b_null

# grouped = groupby(df, [:num_ele])
# all_same = true
# for group in grouped
# 	res1 = group[1, :result]
# 	for row_i in 2:size(group, 1)
# 		res2 = group[row_i, :result]
# 		if !check_similarity(res1, res2)
# 			global all_same = false
# 			println("Results are not the same!")
# 			break
# 		end
# 	end
# end

# all_same
