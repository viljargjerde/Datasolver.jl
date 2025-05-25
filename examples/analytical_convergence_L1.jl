include("basic_setup.jl")
using Datasolver
using JSON
results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")

data_pts_range = [16, 32, 64, 128]
ele_range = 2:6
if isfile(results_file)
	println("Loading existing results from file...")
	result_data = JSON.parsefile(results_file)

	L1_lin_errors = reshape(result_data["L1_lin"], length(ele_range), length(data_pts_range))
	L2_lin_errors = reshape(result_data["L2_lin"], length(ele_range), length(data_pts_range))
	L1_nonlin_errors = reshape(result_data["L1_nonlin"], length(ele_range), length(data_pts_range))
	L2_nonlin_errors = reshape(result_data["L2_nonlin"], length(ele_range), length(data_pts_range))
else
	println("Running solver and storing results...")

	L1_lin_errors = zeros(length(ele_range), length(data_pts_range))
	L2_lin_errors = zeros(length(ele_range), length(data_pts_range))
	L1_nonlin_errors = zeros(length(ele_range), length(data_pts_range))
	L2_nonlin_errors = zeros(length(ele_range), length(data_pts_range))

	for (i_ele, num_ele) in enumerate(ele_range)
		for (j_data, numDataPts) in enumerate(data_pts_range)
			β = 0.01
			n = 1
			bar_length = 1.0
			area = 1.0
			bar_E = 1.0e3

			force_func(α, x) =
				[bar_E * area * ((((1 / 2 * β) * pi^2) * sin((pi * x / bar_length))) *
								 (((3 * α^2) * (((β * pi) * cos((pi * x / bar_length)) / bar_length)^2)) +
								  ((((6 * α) * β) * pi) * cos((pi * x / bar_length)) / bar_length) + 2) / bar_length^2)]

			lin_force = x -> force_func(0.0, x)
			nonlin_force = x -> force_func(1.0, x)

			xs = 0:0.01:bar_length
			strain_limit = maximum(x[1] for x in nonlin_force.(xs)) / bar_E
			dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)

			linear_problem = fixedBarproblem1D(bar_length, area, lin_force, num_ele, 0.0; right_fixed = true)
			nonlinear_problem = fixedBarproblem1D(bar_length, area, nonlin_force, num_ele, 1.0; right_fixed = true)

			lin_result_L1 = Datasolver.NLP_solver(linear_problem, dataset, use_L1_norm = true, worklimit = 100)
			lin_result_L2 = Datasolver.NLP_solver(linear_problem, dataset, use_L1_norm = false, worklimit = 100)
			nonlin_result_L1 = Datasolver.NLP_solver(nonlinear_problem, dataset, use_L1_norm = true, worklimit = 100)
			nonlin_result_L2 = Datasolver.NLP_solver(nonlinear_problem, dataset, use_L1_norm = false, worklimit = 100)

			xs = [p[1] for p in lin_result_L1.Φ]
			u_ans = [β * sin(pi * x / bar_length) for x in xs]

			L1_lin_errors[i_ele, j_data] = calc_reldiff(xs, Datasolver.get_final(lin_result_L1).u, u_ans)
			L2_lin_errors[i_ele, j_data] = calc_reldiff(xs, Datasolver.get_final(lin_result_L2).u, u_ans)
			L1_nonlin_errors[i_ele, j_data] = calc_reldiff(xs, Datasolver.get_final(nonlin_result_L1).u, u_ans)
			L2_nonlin_errors[i_ele, j_data] = calc_reldiff(xs, Datasolver.get_final(nonlin_result_L2).u, u_ans)
		end
	end

	# Save to JSON
	result_data = Dict(
		"L1_lin" => vec(L1_lin_errors),
		"L2_lin" => vec(L2_lin_errors),
		"L1_nonlin" => vec(L1_nonlin_errors),
		"L2_nonlin" => vec(L2_nonlin_errors),
	)

	mkpath(dirname(results_file))
	open(results_file, "w") do f
		JSON.print(f, result_data)
	end
end

# Subtract minimum per cell across methods
function percent_worse_than_best(errors...)
	all_errors = hcat(errors...)
	best_vals = minimum(all_errors, dims = 2)
	return map(e -> ((e .- best_vals) ./ best_vals) * 100, errors)
end

L1_lin_adj, L2_lin_adj, L1_nonlin_adj, L2_nonlin_adj = percent_worse_than_best(
	vec(L1_lin_errors), vec(L2_lin_errors), vec(L1_nonlin_errors), vec(L2_nonlin_errors),
)

reshape_all = x -> reshape(x, size(L1_lin_errors))
L1_lin_adj = reshape_all(L1_lin_adj)
L2_lin_adj = reshape_all(L2_lin_adj)
L1_nonlin_adj = reshape_all(L1_nonlin_adj)
L2_nonlin_adj = reshape_all(L2_nonlin_adj)

p1 = heatmap(data_pts_range, ele_range, L1_lin_adj,
	xlabel = "Data Points", ylabel = "Elements", title = "L1 Linear Error % Worse Than Best")

p2 = heatmap(data_pts_range, ele_range, L2_lin_adj,
	xlabel = "Data Points", ylabel = "Elements", title = "L2 Linear Error % Worse Than Best")

p3 = heatmap(data_pts_range, ele_range, L1_nonlin_adj,
	xlabel = "Data Points", ylabel = "Elements", title = "L1 Nonlinear Error % Worse Than Best")

p4 = heatmap(data_pts_range, ele_range, L2_nonlin_adj,
	xlabel = "Data Points", ylabel = "Elements", title = "L2 Nonlinear Error % Worse Than Best")
plot(p1, p2, p3, p4, layout = (2, 2), size = (800, 600))

plotly()
surface(
	data_pts_range, ele_range, L1_lin_errors;
	xlabel = "Data Points",
	ylabel = "Elements",
	zlabel = "Relative Error",
	title = "Error Surfaces by Method",
	color = :blue,
	label = "L1 Linear",
	legend = :topright,
	colorbar = nothing,
)

surface!(
	data_pts_range, ele_range, L2_lin_errors;
	color = :red,
	label = "L2 Linear",
)

surface!(
	data_pts_range, ele_range, L1_nonlin_errors;
	color = :green,
	label = "L1 Nonlinear",
)

surface!(
	data_pts_range, ele_range, L2_nonlin_errors;
	color = :orange,
	label = "L2 Nonlinear",
)

uncomment_pgfplotsset_blocks(dirname(results_file))
