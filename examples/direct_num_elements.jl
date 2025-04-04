include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables


# num_data_pts = 2^5
num_eles = [2^n for n in 4:12]

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	# Solve a problem to ensure precompilation
	directSolverNonLinearBar(
		linear_problem,
		create_dataset(10, x -> bar_E * x, -strain_limit, strain_limit);
		random_init_data = false,
	)
	for i in 1:100
		for num_ele in num_eles
			linear_problem, _ = get_problems(num_ele)
			for random_init in [false, true]
				t1 = time()
				result = directSolverNonLinearBar(
					linear_problem,
					dataset;
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

process_results(df, results_file)


