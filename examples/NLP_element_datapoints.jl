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
	xlabel = "Elements",
	ylabel = "Data Points",
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

uncomment_pgfplotsset_blocks(dirname(results_file))

if sum(timeout_mask) > 0

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
end
uncomment_pgfplotsset_blocks(dirname(results_file))
