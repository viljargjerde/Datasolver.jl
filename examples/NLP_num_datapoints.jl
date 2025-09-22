include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

numDataPts = [2^n for n in 5:12]
num_ele = 16

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
	for num_data_pts in numDataPts
		dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
		for i in 1:10
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
					"Work" => result.solvetime[1],
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


table = process_results(df, results_file, ("Work", "Work"))
# table = process_results(select(df, Not("Solve time")), results_file, ("Work", "Work"))
plot(scale = :log2, legend = :topleft, xlabel = "Number of datapoints", ylabel = "Solve time [s]")
plot!(table[2:end, "Datapoints"], table[2:end, "Nullspace initialization"], marker = :circle, label = "Nullspace initialization")
plot!(table[2:end, "Datapoints"], table[2:end, "Random initialization"], marker = :circle, label = "Random initialization")
a_null, b_null, f1 = estimate_powerlaw(table[2:end, "Datapoints"], table[2:end, "Nullspace initialization"])
plot!(table[2:end, "Datapoints"], a_null .* table[2:end, "Datapoints"] .^ b_null, label = nothing, linestyle = :dash)
a_rand, b_rand, f2 = estimate_powerlaw(table[2:end, "Datapoints"], table[2:end, "Random initialization"])
plot!(table[2:end, "Datapoints"], a_rand .* table[2:end, "Datapoints"] .^ b_rand, label = nothing, linestyle = :dash)
# annotate!([(256.0, 0.3867922206833348, text("mytext", :red, :right, 3))])



# Define the difference function
diff(x) = f1(x) - f2(x)


# Binary search for the root of diff(x) = 0
function binary_search(f, low, high; tol = 1e-8, max_iter = 1000)
	mid = (low + high) / 2
	for i in 1:max_iter
		mid = (low + high) / 2
		fmid = f(mid)

		if abs(fmid) < tol
			println("Converged after $i iterations")
			return mid
		elseif f(low) * fmid < 0
			high = mid
		else
			low = mid
		end
	end
	error("Binary search did not converge $mid")
end

x_intersect = binary_search(diff, 100, 100000000.0)
y_intersect = f1(x_intersect)

println("Intersection at x = $x_intersect, y = $y_intersect")



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
