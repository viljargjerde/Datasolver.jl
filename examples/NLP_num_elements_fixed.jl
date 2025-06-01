include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

num_eles = [round(Int, 2^n) for n in 3:0.5:6]
push!(num_eles, 2^6 + 1)
# num_eles = [round(Int, 2^n) for n in 3:0.5:6.5]

results_file = joinpath("../master_thesis/both-fixed-figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, _ = get_problems(num_ele, right_fixed = true)
if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")

	for num_ele in num_eles
		linear_problem, _ = get_problems(num_ele, right_fixed = true)
		for random_init in [false, true]
			t1 = time()
			result = NLP_solver(
				linear_problem,
				dataset;
				use_L1_norm = true,
				random_init_data = random_init,
				parameter_file = "NLP_params.prm",
				worklimit = 2000,
			)
			t2 = time()

			push!(results_list, Dict(
				"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
				"Elements" => num_ele,
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


df = DataFrame(Dict.(results_list))

total_work = 0.0
total_time = 0.0
for res in df.Result
	work = res["solvetime"][1]
	solvetime = res["solvetime"][2]
	if work > 0.1 && solvetime > 0.1
		global total_work += work
		global total_time += solvetime
	end
end

println("Work to solvetime ratio: ", total_work / total_time)

table = unstack(select(df, Not(["Result"])), :Initialization, :Work)


xs_fitted = table[3, "Elements"]:table[end, "Elements"]
p = plot(scale = :log2, xlabel = "Number of elements", ylabel = "Work units", palette = paired_colors) # :Paired_12 ,:tableau_20
x_ticks = num_eles[3:2:end]
plot!(table[3:end, "Elements"], table[3:end, "Nullspace initialization"], marker = :circle, label = L"\textbf{Nullspace initialization}", xticks = x_ticks)
a_null, b_null, f1 = estimate_powerlaw(table[3:end-1, "Elements"], table[3:end-1, "Nullspace initialization"])


plot!(xs_fitted, f1.(xs_fitted), label = L"W(\mathtt{M}) = %$(latex_sci(a_null)) \ \mathtt{M}^{%$(round(b_null, sigdigits = 2))}", linestyle = :dash)
a_rand, b_rand, f2 = estimate_powerlaw(table[3:end-1, "Elements"], table[3:end-1, "Random initialization"])


plot!(table[3:end, "Elements"], table[3:end, "Random initialization"], marker = :circle, label = L"\textbf{No initialization}")
plot!(xs_fitted, f2.(xs_fitted), label = L"W(\mathtt{M}) = %$(latex_sci(a_null)) \ \mathtt{M}^{%$(round(b_null, sigdigits = 2))}", linestyle = :dash)

savefig(replace(results_file, "results.json" => "lineplot.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))

update_tex_command(all_results_file, "NLPFixedElementsPowerlawNull", format(FormatExpr("W(\\mathtt{{M}}) \\propto \\mathtt{{M}}^{{{:.2f}}}"), b_null))
update_tex_command(all_results_file, "NLPFixedElementsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))
update_tex_command(all_results_file, "NLPFixedElementsPowerlawRand", format(FormatExpr("W(\\mathtt{{M}}) \\propto \\mathtt{{M}}^{{{:.2f}}}"), b_rand))
update_tex_command(all_results_file, "NLPFixedElementsSpeedRatio", format(FormatExpr("{:.1f}"), a_rand / a_null))

p
