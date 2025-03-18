include("../../basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

numDataPts = [2^n for n in 5:10]
num_ele = 12

results_file = "examples/NLP_solver_comparisons/initialization_num_datapoints/results.json"
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
				)
				t2 = time()

				push!(results_list, Dict(
					"random_init" => random_init,
					"num_data_pts" => num_data_pts,
					"solve_time" => t2 - t1,
					"result" => result,
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

# Compute statistics
grouped = combine(groupby(df, [:num_data_pts, :random_init]),
	:solve_time => median => :median_solve_time,
	:solve_time => (x -> quantile(x, 0.25)) => :q25,
	:solve_time => (x -> quantile(x, 0.75)) => :q75,
)


table = unstack(grouped, :num_data_pts, :random_init, :median_solve_time)
rename!(table, Dict(:num_data_pts => "Data points"))
rename!(table, Dict("true" => "Random"))
rename!(table, Dict("false" => "Non-Random"))

select!(table, "Data points", "Random", "Non-Random")
open(replace(results_file, ".json" => ".tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
	pretty_table(table, show_subheader = false)
end

# Extract error bars for plotting
error_random = grouped.q75[grouped.random_init.==true] .- grouped.q25[grouped.random_init.==true]
error_nonrandom = grouped.q75[grouped.random_init.==false] .- grouped.q25[grouped.random_init.==false]

table_random = table."Random"
table_nonrandom = table."Non-Random"
data_points = table."Data points"

plot(data_points, table_random, label = "Random initialization", marker = :circle, scale = :log2, yerror = error_random)
plot!(data_points, table_nonrandom, label = "Non-Random initialization", marker = :square, ylabel = "Median Solve Time (s)", xlabel = "Number of data points", yerror = error_nonrandom)

savefig(replace(results_file, ".json" => ".png"))

grouped = groupby(df, [:num_data_pts])
all_same = true
for group in grouped
	res1 = group[1, :result]
	for row_i in 2:size(group, 1)
		res2 = group[row_i, :result]
		if !check_similarity(res1, res2)
			global all_same = false
			println("Results are not the same!")
			break
		end
	end
end

all_same
