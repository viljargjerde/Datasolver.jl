include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using ColorSchemes

numDataPts = [2^n + 1 for n in 4:10]

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")

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
				parameter_file = "NLP_params.prm")
			t2 = time()

			push!(results_list, Dict(
				"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
				"Datapoints" => num_data_pts,
				"Solve time" => t2 - t1,
				"Work" => result.solvetime[1],
				"Result" => result,
			))
			# Save results to file
			open(results_file, "w") do f
				JSON.print(f, results_list)
			end
		end
	end
end
# # Convert results to DataFrame
df = DataFrame(Dict.(results_list))

# table = process_results(df, results_file, ("Work", "Work"))
table = unstack(select(df, Not(["Solve time", "Result"])), :Initialization, :Work)
# table = process_results(select(df, Not("Solve time")), results_file, ("Work", "Work"))
p = plot(scale = :log2, xlabel = "Data points", ylabel = "Work units", palette = paired_colors) # :Paired_12 ,:tableau_20
plot!(table[2:end, "Datapoints"], table[2:end, "Nullspace initialization"], marker = :circle, label = L"\textbf{Nullspace initialization}")
a_null, b_null, f2 = estimate_powerlaw(table[2:end, "Datapoints"], table[2:end, "Nullspace initialization"])
update_tex_command(all_results_file, "NLPDatapointsPowerlawNull", format(FormatExpr("W(\\mathtt{{D}}) \\propto \\mathtt{{D}}^{{{:.2f}}}"), b_null))
# update_tex_command(all_results_file, "NLPDatapointsPowerlawNull", format(FormatExpr("y = {:.2g}D^{{{:.2f}}}"), a_null, b_null))
update_tex_command(all_results_file, "NLPDatapointsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))

plot!(table[2:end, "Datapoints"], a_null .* table[2:end, "Datapoints"] .^ b_null, label = L"W(\mathtt{D}) = %$(latex_sci(a_null)) \ \mathtt{D}^{%$(round(b_null, sigdigits = 2))}", linestyle = :dash)
a_rand, b_rand, f2 = estimate_powerlaw(table[2:end, "Datapoints"], table[2:end, "Random initialization"])
update_tex_command(all_results_file, "NLPDatapointsPowerlawRand", format(FormatExpr("W(\\mathtt{{D}}) \\propto \\mathtt{{D}}^{{{:.2f}}}"), b_rand))
update_tex_command(all_results_file, "NLPDatapointsSpeedRatio", format(FormatExpr("{:.1f}"), a_rand / a_null))
# update_tex_command(all_results_file, "NLPDatapointsPowerlawRand", format(FormatExpr("y = {:.2f}x^{{{:.2f}}}"), a_rand, b_rand))

plot!(table[2:end, "Datapoints"], table[2:end, "Random initialization"], marker = :circle, label = L"\textbf{No initialization}")
plot!(table[2:end, "Datapoints"], a_rand .* table[2:end, "Datapoints"] .^ b_rand, label = L"W(\mathtt{D}) = %$(latex_sci(a_rand)) \ \mathtt{D}^{%$(round(b_rand, sigdigits = 2))}", linestyle = :dash)
savefig(replace(results_file, "results.json" => "lineplot.tex"))

uncomment_pgfplotsset_blocks(dirname(results_file))

p
