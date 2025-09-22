include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables


# num_data_pts = 2^5
num_eles = [2^n for n in 6:11]
dataset = create_dataset(2^12, x -> bar_E * x, 0.0, 2 * strain_limit)

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
	for num_ele in num_eles
		@show num_ele
		for i in 1:100
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

table = process_results(df, results_file)

xs_fitted = table[3, "Elements"]:table[end, "Elements"]
p = plot(scale = :log2, xlabel = "Number of elements", ylabel = "Mean solve time (s)", palette = paired_colors) # :Paired_12 ,:tableau_20
x_ticks = num_eles[3:2:end]

plot!(table[3:end, "Elements"], table[3:end, "Nullspace initialization"], marker = :circle, label = L"\textbf{Nullspace initialization}", xticks = x_ticks)
a_null, b_null, f1 = estimate_powerlaw(table[3:end-1, "Elements"], table[3:end-1, "Nullspace initialization"])
plot!(xs_fitted, f1.(xs_fitted), label = L"W(\mathtt{M}) = %$(latex_sci(a_null)) \ \mathtt{M}^{%$(round(b_null, sigdigits = 2))}", linestyle = :dash)


a_rand, b_rand, f2 = estimate_powerlaw(table[3:end-1, "Elements"], table[3:end-1, "Random initialization"])
plot!(table[3:end, "Elements"], table[3:end, "Random initialization"], marker = :circle, label = L"\textbf{Random initialization}")
plot!(xs_fitted, f2.(xs_fitted), label = L"W(\mathtt{M}) = %$(latex_sci(a_rand)) \ \mathtt{M}^{%$(round(b_rand, sigdigits = 2))}", linestyle = :dash)

savefig(replace(results_file, "results.json" => "lineplot.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))


update_tex_command(all_results_file, "DirectElementsPowerlawNull", format(FormatExpr("T(\\mathtt{{M}}) \\propto \\mathtt{{M}}^{{{:.2f}}}"), b_null))
update_tex_command(all_results_file, "DirectElementsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))

update_tex_command(all_results_file, "DirectElementsSpeedRatio", format(FormatExpr("{:.1f}"), a_rand / a_null))
update_tex_command(all_results_file, "DirectElementsPowerlawRand", format(FormatExpr("T(\\mathtt{{M}}) \\propto \\mathtt{{M}}^{{{:.2f}}}"), b_rand))

p
