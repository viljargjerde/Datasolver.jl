include("basic_setup.jl")
using Datasolver

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

# nonlin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-8)
# plot_results(lin_result, dataset = dataset, title = "Linear solver")
# plot_results(nonlin_result, dataset = dataset, title = "Nonlinear solver")
# # result = Datasolver.greedyLocalSearchSolverNonLinearBar(linear_problem, dataset)
# @show lin_result.cost[end]
# @show nonlin_result.cost[end]
# # result = NLP_solver(linear_problem, dataset, use_L1_norm = false, worklimit = 100)

# @show L2_nonlin = round(calc_reldiff(xs, Datasolver.get_final(nonlin_result).u, u_ans), digits = 4)


# plot(xs, u_ans, label = "Analytical solution", color = :black, linewidth = 2)
# plot!(xs, get_final(lin_result_L1).u, label = "L1", color = :blue, linewidth = 2)
# plot!(xs, get_final(lin_result_L2).u, label = "L2", color = :red, linewidth = 2)

# ################# Linear #####################
# N_datapoints = [2^n for n in 2:10]
# N_elements = [2^n for n in 2:10]
# analytical_u = []
# results = Vector{NamedTuple}()
# for N_d in N_datapoints, N_e in N_elements
# 	@show N_d, N_e
# 	local dataset = create_dataset(N_d, x -> bar_E * x, -strain_limit, strain_limit)
# 	local linear_problem = fixedBarproblem1D(
# 		bar_length,
# 		area,
# 		lin_force,
# 		N_e,
# 		0.0;
# 		right_fixed = true,
# 	)
# 	local result = greedyLocalSearchSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-9)
# 	push!(results, (N_datapoints = N_d, N_elements = N_e, result = result))

# 	xs = [p[1] for p in result.Φ]
# 	u_ans = [β * sin(pi * x / bar_length) for x in xs]
# 	push!(analytical_u, u_ans)
# end
# Datasolver.convergence_analysis(results, analytical_u)


# ################# Nonlinear #####################

# nonlin_analytical_u = []
# nonlin_results = Vector{NamedTuple}()
# for N_d in N_datapoints, N_e in N_elements
# 	@show N_d, N_e
# 	local dataset = create_dataset(N_d, x -> bar_E * x, -strain_limit, strain_limit)
# 	local nonlinear_problem = fixedBarproblem1D(
# 		bar_length,
# 		area,
# 		nonlin_force,
# 		N_e,
# 		1.0;
# 		right_fixed = true,
# 	)
# 	local result = greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-9)
# 	push!(nonlin_results, (N_datapoints = N_d, N_elements = N_e, result = result))

# 	xs = [p[1] for p in result.Φ]
# 	u_ans = [β * sin(pi * x / bar_length) for x in xs]
# 	push!(nonlin_analytical_u, u_ans)
# end
# Datasolver.convergence_analysis(nonlin_results, nonlin_analytical_u)


# Datasolver.convergence_analysis_dual(results, analytical_u, nonlin_results, nonlin_analytical_u)

# savefig(replace(results_file, "results.json" => "convergence_contour.tex"))

# plot_results(result, dataset = dataset)

# plot_results(result, dataset = dataset)


# gr()
# ########################################################################################################
# # Right end free
# ########################################################################################################
# begin
# 	numDataPts = 2048
# 	num_ele = 16
# 	f_0 = 200.0
# 	# Force function
# 	force = x -> [f_0]

# 	# Maximum strain limit for dataset generation
# 	strain_limit = f_0 * bar_length / (bar_E * area)

# 	# Dataset creation
# 	dataset = create_dataset(numDataPts, x -> bar_E * x, 0.0, 2 * strain_limit)

# 	# Problem definitions
# 	linear_problem = fixedBarproblem1D(
# 		bar_length, area, force, num_ele, 0.0; right_fixed = false,
# 	)

# 	nonlinear_problem = fixedBarproblem1D(
# 		bar_length, area, force, num_ele, 1.0; right_fixed = false,
# 	)

# 	# Solvers
# 	# lin_result = Datasolver.NLP_solver(linear_problem, dataset, use_L1_norm = false, worklimit = 1000)
# 	# nonlin_result = Datasolver.NLP_solver(nonlinear_problem, dataset, use_L1_norm = false, worklimit = 1000, random_init_data = true)


# 	lin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-8)
# 	nonlin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-8)
# 	# Visualization
# 	plot_results(lin_result, dataset = dataset, title = "Linear solver")
# 	plot_results(nonlin_result, dataset = dataset, title = "Nonlinear solver")

# 	# Cost
# 	@show lin_result.cost[end]
# 	@show nonlin_result.cost[end]

# 	# Extract positions
# 	xs = [p[1] for p in lin_result.Φ]

# 	# Analytical solutions
# 	u_lin(x) = -f_0 / (2 * bar_E * area) * x^2 + f_0 * bar_length / (bar_E * area) * x
# 	# u_nonlin(x) = -(bar_E / f_0) * (sqrt(1 + 2 * f_0 * (bar_length - x) / bar_E) + sqrt(1 + 2 * f_0 / bar_E))

# 	function u_nonlin(x)
# 		coeff = (2 * f_0) / (bar_E * area)
# 		term1 = (1 + coeff * bar_length)^(3 / 2)
# 		term2 = (1 + coeff * (bar_length - x))^(3 / 2)
# 		return (bar_E / (3 * f_0 / area)) * (term1 - term2) - x
# 	end



# 	xs_analytical = range(xs[1], stop = xs[end], length = 100)

# 	# Relative errors
# 	@show L2_lin = round(calc_reldiff(xs, Datasolver.get_final(lin_result).u, u_lin.(xs)), digits = 4)
# 	@show L2_nonlin = round(calc_reldiff(xs, Datasolver.get_final(nonlin_result).u, u_nonlin.(xs)), digits = 4)

# 	# Plot comparison
# 	plot(palette = paired_colors)
# 	plot!(xs_analytical, u_lin.(xs_analytical), label = "Analytical Linear", linewidth = 2)
# 	plot!(xs, get_final(lin_result).u, label = "Data-driven Linear", linestyle = :dash, linewidth = 2)
# 	plot!(xs_analytical, u_nonlin.(xs_analytical), label = "Analytical NonLinear", linewidth = 2)
# 	plot!(xs, get_final(nonlin_result).u, label = "Data-driven NonLinear", linestyle = :dash, linewidth = 2)



# end



# # @code_warntype Datasolver.directSolverNonLinearBar(linear_problem, dataset, init)
