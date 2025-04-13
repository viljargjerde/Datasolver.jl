include("basic_setup.jl")
using Datasolver


gr()

numDataPts = 1024
num_ele = 128
f_0 = 100.0
n = 1
force = x -> [f_0 * sin(n * pi * x / bar_length)]

strain_limit = f_0 * linear_problem.length / (bar_E * linear_problem.area);

dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)

linear_problem = fixedBarproblem1D(
	bar_length,
	area,
	force,
	num_ele,
	0.0;
	right_fixed = true,
)

nonlinear_problem = fixedBarproblem1D(
	bar_length,
	area,
	force,
	num_ele,
	1.0;
	right_fixed = true,
)
lin_result = Datasolver.directSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-8)
nonlin_result = Datasolver.directSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-8)
plot_results(lin_result, dataset = dataset, title = "Linear solver")
plot_results(nonlin_result, dataset = dataset, title = "Nonlinear solver")
# result = Datasolver.greedyLocalSearchSolverNonLinearBar(linear_problem, dataset)
@show lin_result.cost[end]
@show nonlin_result.cost[end]
# result = NLP_solver(linear_problem, dataset, use_L1_norm = false, worklimit = 100)

xs = [p[1] for p in lin_result.Φ]
α = f_0 / (bar_E * area * (n * pi / bar_length)^2)
u_ans = [α * sin(n * pi * x / bar_length) for x in xs]
@show L2_lin = round(calc_reldiff(xs, Datasolver.get_final(lin_result).u, u_ans), digits = 4)
@show L2_nonlin = round(calc_reldiff(xs, Datasolver.get_final(nonlin_result).u, u_ans), digits = 4)


plot(xs, u_ans, label = "Analytical solution", color = :black, linewidth = 2)
plot!(xs, get_final(lin_result).u, label = "Linear solver", color = :blue, linewidth = 2)
plot!(xs, get_final(nonlin_result).u, label = "NonLinear solver", color = :red, linewidth = 2)


N_datapoints = [2^n for n in 5:11]
N_elements = [2^n for n in 5:11]
analytical_u = []
results = Vector{NamedTuple}()
for N_d in N_datapoints, N_e in N_elements
	@show N_d, N_e
	local dataset = create_dataset(N_d, x -> bar_E * x, -strain_limit, strain_limit)
	local linear_problem = fixedBarproblem1D(
		bar_length,
		area,
		force,
		N_e,
		0.0;
		right_fixed = true,
	)
	local result = directSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-9)
	push!(results, (N_datapoints = N_d, N_elements = N_e, result = result))

	xs = [p[1] for p in result.Φ]
	u_ans = [α * sin(n * pi * x / bar_length) for x in xs]
	push!(analytical_u, u_ans)
end
Datasolver.convergence_analysis(results, analytical_u)

# plot_results(result, dataset = dataset)

plot_results(result, dataset = dataset)


gr()
########################################################################################################
# Right end free
########################################################################################################
begin
	numDataPts = 2048
	num_ele = 16
	f_0 = 200.0
	# Force function
	force = x -> [f_0]

	# Maximum strain limit for dataset generation
	strain_limit = f_0 * bar_length / (bar_E * area)

	# Dataset creation
	dataset = create_dataset(numDataPts, x -> bar_E * x, 0.0, 2 * strain_limit)

	# Problem definitions
	linear_problem = fixedBarproblem1D(
		bar_length, area, force, num_ele, 0.0; right_fixed = false,
	)

	nonlinear_problem = fixedBarproblem1D(
		bar_length, area, force, num_ele, 1 / 3; right_fixed = false,
	)

	# Solvers
	# lin_result = Datasolver.NLP_solver(linear_problem, dataset, use_L1_norm = false, worklimit = 1000)
	# nonlin_result = Datasolver.NLP_solver(nonlinear_problem, dataset, use_L1_norm = false, worklimit = 1000, random_init_data = true)


	lin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-8)
	nonlin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-8)
	# Visualization
	plot_results(lin_result, dataset = dataset, title = "Linear solver")
	plot_results(nonlin_result, dataset = dataset, title = "Nonlinear solver")

	# Cost
	@show lin_result.cost[end]
	@show nonlin_result.cost[end]

	# Extract positions
	xs = [p[1] for p in lin_result.Φ]

	# Analytical solutions
	u_lin(x) = -f_0 / (2 * bar_E * area) * x^2 + f_0 * bar_length / (bar_E * area) * x
	# u_nonlin(x) = -(bar_E / f_0) * (sqrt(1 + 2 * f_0 * (bar_length - x) / bar_E) + sqrt(1 + 2 * f_0 / bar_E))

	function u_nonlin(x)
		coeff = (2 * f_0) / (bar_E * area)
		term1 = (1 + coeff * bar_length)^(3 / 2)
		term2 = (1 + coeff * (bar_length - x))^(3 / 2)
		return (bar_E / (3 * f_0 / area)) * (term1 - term2) - x
	end



	xs_analytical = range(xs[1], stop = xs[end], length = 100)

	# Relative errors
	@show L2_lin = round(calc_reldiff(xs, Datasolver.get_final(lin_result).u, u_lin.(xs)), digits = 4)
	@show L2_nonlin = round(calc_reldiff(xs, Datasolver.get_final(nonlin_result).u, u_nonlin.(xs)), digits = 4)

	# Plot comparison
	plot(palette = paired_colors)
	plot!(xs_analytical, u_lin.(xs_analytical), label = "Analytical Linear", linewidth = 2)
	plot!(xs, get_final(lin_result).u, label = "Data-driven Linear", linestyle = :dash, linewidth = 2)
	plot!(xs_analytical, u_nonlin.(xs_analytical), label = "Analytical NonLinear", linewidth = 2)
	plot!(xs, get_final(nonlin_result).u, label = "Data-driven NonLinear", linestyle = :dash, linewidth = 2)



end



# @code_warntype Datasolver.directSolverNonLinearBar(linear_problem, dataset, init)
