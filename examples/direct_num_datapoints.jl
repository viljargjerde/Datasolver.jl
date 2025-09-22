include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

numDataPts = [2^n for n in 5:16]
num_ele = 256
results_file = results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
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
		for num_data_pts in numDataPts
			dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
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
using LaTeXStrings


table = process_results(df, results_file, ("Solve time", "Median solve time (s)"))
p = plot(scale = :log2, xlabel = "Data points", ylabel = "Median solve time (s)", palette = paired_colors) # :Paired_12 ,:tableau_20
plot!(table[1:end, "Datapoints"], table[1:end, "Nullspace initialization"], marker = :circle, label = L"\textbf{Nullspace initialization}")
a_null, b_null, f2 = estimate_powerlaw(table[1:end, "Datapoints"], table[1:end, "Nullspace initialization"])
plot!(table[1:end, "Datapoints"], a_null .* table[1:end, "Datapoints"] .^ b_null, label = L"W(D) = %$(latex_sci(a_null)) \ D^{%$(round(b_null, sigdigits = 2))}", linestyle = :dash)



a_rand, b_rand, f2 = estimate_powerlaw(table[1:end, "Datapoints"], table[1:end, "Random initialization"])

plot!(table[1:end, "Datapoints"], table[1:end, "Random initialization"], marker = :circle, label = L"\textbf{Random initialization}") # "Random initialization")
plot!(table[1:end, "Datapoints"], a_rand .* table[1:end, "Datapoints"] .^ b_rand, label = L"W(D) = %$(latex_sci(a_rand)) D^{%$(round(b_rand, sigdigits = 2))}", linestyle = :dash)

savefig(replace(results_file, "results.json" => "lineplot.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))


update_tex_command(all_results_file, "DirectDatapointsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))
update_tex_command(all_results_file, "DirectDatapointsPowerlawRand", format(FormatExpr("T(D) \\propto D^{{{:.2f}}}"), b_rand))
update_tex_command(all_results_file, "DirectDatapointsPowerlawNull", format(FormatExpr("T(D) \\propto D^{{{:.2f}}}"), b_null))

p
