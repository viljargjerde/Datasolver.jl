include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

num_data_pts = 2^5
num_eles = [2^n for n in 2:6]

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
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
	for i in 1:100
		for num_ele in num_eles
			linear_problem, _ = get_problems(num_ele)
			dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
			for random_init in [false, true]
				t1 = time()
				result = NLP_solver(
					linear_problem,
					dataset;
					use_L1_norm = true,
					random_init_data = random_init,
				)
				t2 = time()

				push!(results_list, Dict(
					"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
					"Elements" => num_ele,
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
df = DataFrame(results_list)

table = process_results(df, results_file)
a_null, b_null, f1 = estimate_powerlaw(table[3:end, "Elements"], table[3:end, "Nullspace initialization"])
a_rand, b_rand, f2 = estimate_powerlaw(table[3:end, "Elements"], table[3:end, "Random initialization"])
f1(12)
f1(200)
f2(12)
f2(200)
b_rand
b_null

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
