include("basic_setup.jl")
using Plots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables

num_eles = [2^n for n in 2:6]

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []
linear_problem, _ = get_problems(num_ele)
if isfile(results_file)
    println("Loading existing results from file...")
    results_list = JSON.parsefile(results_file)
else
    println("Running solver and storing results...")

    for num_ele in num_eles
        linear_problem, _ = get_problems(num_ele)
        for random_init in [false, true]
            t1 = time()
            result = NLP_solver(
                linear_problem,
                dataset;
                use_L1_norm=true,
                random_init_data=random_init,
            )
            t2 = time()

            push!(results_list, Dict(
                "Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
                "Elements" => num_ele,
                "Work" => result.solvetime[1],
                "Result" => result,
            ))
        end
    end
    # Save results to file
    open(results_file, "w") do f
        JSON.print(f, results_list)
    end
end

df = DataFrame(Dict.(results_list))



# table = process_results(df, results_file, ("Work", "Work"))
table = unstack(select(df, Not(["Result"])), :Initialization, :Work)
# table = process_results(select(df, Not("Solve time")), results_file, ("Work", "Work"))
p = plot(scale=:log2, legend=:topleft, xlabel="Number of elements", ylabel="Work units", palette=paired_colors)
plot!(table[2:end, "Elements"], table[2:end, "Nullspace initialization"], marker=:circle, label="Nullspace initialization")
a_null, b_null, f2 = estimate_powerlaw(table[2:end, "Elements"], table[2:end, "Nullspace initialization"])
update_tex_command(all_results_file, "NLPElementsPowerlawNull", format(FormatExpr("y \\propto x^{{{:.2f}}}"), b_null))
# update_tex_command(all_results_file, "NLPElementsPowerlawNull", format(FormatExpr("y = {:.2g}x^{{{:.2f}}}"), a_null, b_null))
update_tex_command(all_results_file, "NLPElementsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))

plot!(table[2:end, "Elements"], a_null .* table[2:end, "Elements"] .^ b_null, label=nothing, linestyle=:dash)
a_rand, b_rand, f2 = estimate_powerlaw(table[2:end, "Elements"], table[2:end, "Random initialization"])
update_tex_command(all_results_file, "NLPElementsPowerlawRand", format(FormatExpr("y \\propto x^{{{:.2f}}}"), b_rand))
# update_tex_command(all_results_file, "NLPElementsPowerlawRand", format(FormatExpr("y = {:.2f}x^{{{:.2f}}}"), a_rand, b_rand))

plot!(table[2:end, "Elements"], table[2:end, "Random initialization"], marker=:circle, label="Random initialization")
plot!(table[2:end, "Elements"], a_rand .* table[2:end, "Elements"] .^ b_rand, label=nothing, linestyle=:dash)
savefig(replace(results_file, "results.json" => "lineplot.tex"))
p
# # Convert results to DataFrame
# df = DataFrame(results_list)

# table = process_results(df, results_file)
# a_null, b_null, f1 = estimate_powerlaw(table[3:end, "Elements"], table[3:end, "Nullspace initialization"])
# a_rand, b_rand, f2 = estimate_powerlaw(table[3:end, "Elements"], table[3:end, "Random initialization"])
# f1(12)
# f1(200)
# f2(12)
# f2(200)
# b_rand
# b_null

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
