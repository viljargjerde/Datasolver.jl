include("basic_setup.jl")

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")

linear_problem = fixedBarproblem1D(
	bar_length,
	area,
	force,
	num_ele,
	0.0;
	right_fixed = true,
)
dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
s = Datasolver.get_initialization_s(linear_problem)

dataset.S[4] = minimum(s)

best_idxs = Datasolver.find_closest_idx(dataset.S, s)
S = dataset.S[best_idxs]
E = dataset.E[best_idxs]




if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	results_list = []
	push!(results_list, NLP_solver(
		linear_problem,
		dataset;
		use_L1_norm = false,
		random_init_data = false,
		parameter_file = "NLP_params.prm",
		worklimit = 100,
	))
	push!(results_list, directSolverNonLinearBar(
		linear_problem,
		dataset;
		random_init_data = false,
	))
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end




NLP_results = results_list[1]
direct_results = results_list[2]


scatter(dataset.E, dataset.S, label = "Data", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, line = :dash, markercolor = paired_colors[1])

scatter!(E, S, label = "Initialization", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :square, markercolor = paired_colors[5], line = :dash)


scatter!(direct_results["e"][end], direct_results["s"][end], label = "ADM", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :circle, markercolor = paired_colors[12], line = :dash)


scatter!(NLP_results["e"][end], NLP_results["s"][end], label = "MINLP", xlabel = "Strain", ylabel = "Stress", legend = :topleft, marker = :star, markersize = 2, markercolor = paired_colors[9], line = :dash)

# savefig(replace(results_file, "results.json" => "example-solution.tex"));
uncomment_pgfplotsset_blocks(dirname(results_file))
