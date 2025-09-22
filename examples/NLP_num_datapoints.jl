include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

numDataPts = [2^n for n in 5:10]
num_ele = 12

joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	# Solve a problem to ensure precompilation
	NLP_solver(
		linear_problem,
		create_dataset(10, x -> bar_E * x, -strain_limit, strain_limit);
		use_L1_norm = true,
		random_init_data = false,
	)
	for i in 1:10
		for num_data_pts in numDataPts
			dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
			for random_init in [false, true]
				t1 = time()
				result = NLP_solver(
					linear_problem,
					dataset;
					use_L1_norm = true,
					random_init_data = random_init,
					verbose = false,
				)
				t2 = time()

				push!(results_list, Dict(
					"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
					"Datapoints" => num_data_pts,
					"Solve time" => t2 - t1,
					"Result" => result,
				))
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


process_results(df, results_file)
# grouped = groupby(df, ["Datapoints"])
# all_same = true
# for group in grouped
# 	# @show group
# 	res1 = group[1, :Result]
# 	for row_i in 2:size(group, 1)
# 		res2 = group[row_i, :Result]
# 		if !check_similarity(res1, res2)
# 			global all_same = false
# 			break
# 		end
# 	end
# end

# all_same
