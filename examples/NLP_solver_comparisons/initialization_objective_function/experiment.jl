include("../../basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

num_ele = 20
results_file = "examples/NLP_solver_comparisons/initialization_objective_function/results.json"
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
		dataset;
		use_L1_norm = true,
		random_init_data = false,
	)
	linear_problem, nonlinear_problem = get_problems(num_ele)
	for i in 1:10
		for is_L1 in [true, false]

			for random_init in [false, true]
				t1 = time()
				result = NLP_solver(
					linear_problem,
					dataset;
					use_L1_norm = is_L1,
					random_init_data = random_init,
				)
				t2 = time()

				push!(results_list, Dict(
					"random_init" => random_init,
					"is_L1" => is_L1,
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
df = DataFrame(results_list)

# Compute statistics
grouped = combine(groupby(df, [:is_L1, :random_init]),
	:solve_time => median => :median_solve_time,
	:solve_time => (x -> quantile(x, 0.25)) => :q25,
	:solve_time => (x -> quantile(x, 0.75)) => :q75,
)


table = unstack(grouped, :is_L1, :random_init, :median_solve_time)
table.is_L1 = replace(table.is_L1, false => "L2", true => "L1")
rename!(table, Dict(:is_L1 => "Obective function"))
rename!(table, Dict("true" => "Random"))
rename!(table, Dict("false" => "Non-Random"))

select!(table, "Obective function", "Random", "Non-Random")
open(replace(results_file, ".json" => ".tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
	pretty_table(table, show_subheader = false)
end

# Extract error bars for plotting
error_random = grouped.q75[grouped.random_init.==true] .- grouped.q25[grouped.random_init.==true]
error_nonrandom = grouped.q75[grouped.random_init.==false] .- grouped.q25[grouped.random_init.==false]

table_random = table."Random"
table_nonrandom = table."Non-Random"
strain_measure = table."Obective function"

plot(strain_measure, table_random, label = "Random initialization", marker = :circle, scale = :log2, yerror = error_random)
plot!(strain_measure, table_nonrandom, label = "Non-Random initialization", marker = :square, ylabel = "Median Solve Time (s)", xlabel = "Number of Elements", yerror = error_nonrandom)

savefig(replace(results_file, ".json" => ".png"))


