# include("basic_setup.jl")
# using Plots
# using DataFrames
# using Statistics
# using CSV
# using JSON
# using PrettyTables
# using ColorSchemes

# numDataPts = [2^n + 1 for n in 5:10]

# results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
# results_list = []
# if isfile(results_file)
# 	println("Loading existing results from file...")
# 	results_list = JSON.parsefile(results_file)
# else
# 	println("Running solver and storing results...")

# 	for num_data_pts in numDataPts
# 		local dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit)
# 		for random_init in [false, true]
# 			t1 = time()
# 			result = NLP_solver(
# 				nonlinear_problem,
# 				dataset;
# 				use_L1_norm = true,
# 				random_init_data = random_init,
# 				verbose = false,
# 				parameter_file = "NLP_params.prm")
# 			t2 = time()

# 			push!(results_list, Dict(
# 				"Initialization" => random_init ? "Random initialization" : "Nullspace initialization",
# 				"Datapoints" => num_data_pts,
# 				"Solve time" => t2 - t1,
# 				"Work" => result.solvetime[1],
# 				"Result" => result,
# 			))
# 		end
# 	end
# 	# Save results to file
# 	open(results_file, "w") do f
# 		JSON.print(f, results_list)
# 	end
# end
# # # Convert results to DataFrame
# df = DataFrame(Dict.(results_list))



# # table = process_results(df, results_file, ("Work", "Work"))
# table = unstack(select(df, Not(["Solve time", "Result"])), :Initialization, :Work)
# # table = process_results(select(df, Not("Solve time")), results_file, ("Work", "Work"))
# p = plot(scale = :log2, legend = :topleft, xlabel = "Data points", ylabel = "Work units", palette = paired_colors) # :Paired_12 ,:tableau_20
# plot!(table[2:end, "Datapoints"], table[2:end, "Nullspace initialization"], marker = :circle, label = "Nullspace initialization")
# a_null, b_null, f2 = estimate_powerlaw(table[2:end, "Datapoints"], table[2:end, "Nullspace initialization"])
# update_tex_command(all_results_file, "NLPDatapointsPowerlawNull", format(FormatExpr("t(D) \\propto D^{{{:.2f}}}"), b_null))
# # update_tex_command(all_results_file, "NLPDatapointsPowerlawNull", format(FormatExpr("y = {:.2g}x^{{{:.2f}}}"), a_null, b_null))
# update_tex_command(all_results_file, "NLPDatapointsPowerlawNullB", format(FormatExpr("{:.2g}"), b_null))

# plot!(table[2:end, "Datapoints"], a_null .* table[2:end, "Datapoints"] .^ b_null, label = nothing, linestyle = :dash)
# a_rand, b_rand, f2 = estimate_powerlaw(table[2:end, "Datapoints"], table[2:end, "Random initialization"])
# update_tex_command(all_results_file, "NLPDatapointsPowerlawRand", format(FormatExpr("t(D) \\propto D^{{{:.2f}}}"), b_rand))
# # update_tex_command(all_results_file, "NLPDatapointsPowerlawRand", format(FormatExpr("y = {:.2f}x^{{{:.2f}}}"), a_rand, b_rand))

# plot!(table[2:end, "Datapoints"], table[2:end, "Random initialization"], marker = :circle, label = "No initialization")
# plot!(table[1:end, "Datapoints"], a_rand .* table[1:end, "Datapoints"] .^ b_rand, label = nothing, linestyle = :dash)
# savefig(replace(results_file, "results.json" => "lineplot.tex"))
# p
# # # annotate!([(256.0, 0.3867922206833348, text("mytext", :red, :right, 3))])
# uncomment_pgfplotsset_blocks(dirname(results_file))



# # # Define the difference function
# # diff(x) = f1(x) - f2(x)


# # # Binary search for the root of diff(x) = 0
# # function binary_search(f, low, high; tol = 1e-8, max_iter = 1000)
# # 	mid = (low + high) / 2
# # 	for i in 1:max_iter
# # 		mid = (low + high) / 2
# # 		fmid = f(mid)

# # 		if abs(fmid) < tol
# # 			println("Converged after $i iterations")
# # 			return mid
# # 		elseif f(low) * fmid < 0
# # 			high = mid
# # 		else
# # 			low = mid
# # 		end
# # 	end
# # 	error("Binary search did not converge $mid")
# # end

# # x_intersect = binary_search(diff, 100, 100000000.0)
# # y_intersect = f1(x_intersect)

# # println("Intersection at x = $x_intersect, y = $y_intersect")



# # grouped = groupby(df, ["Datapoints"])
# # all_same = true
# # for group in grouped
# # 	# @show group
# # 	res1 = group[1, :Result]
# # 	for row_i in 2:size(group, 1)
# # 		res2 = group[row_i, :Result]
# # 		if !check_similarity(res1, res2)
# # 			global all_same = false
# # 			break
# # 		end
# # 	end
# # end

# # all_same
