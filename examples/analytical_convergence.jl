include("basic_setup.jl")
using Datasolver

results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")


numDataPts = 32
num_ele = 4
β = 0.01
# f_0 = 100.0
n = 1
bar_length = 1.0
area = 1.0
bar_E = 1.0e3
# force = x -> [f_0 * sin(n * pi * x / bar_length)]
force_func(α, x) =
	[bar_E * area * ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_length))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_length)) / bar_length))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_length)) / bar_length) + 2) / bar_length^(2))]
# force_func(α, x) = [(α / 2 + 1) * (2 * α * π / bar_length * β * cos(π * x / bar_length) + 1) * π^2 / bar_length^2 * β * sin(π * x / bar_length)]
lin_force = x -> force_func(0.0, x)
nonlin_force = x -> force_func(1.0, x)
xs = 0:0.01:bar_length
# plot(xs, [x[1] for x in lin_force.(xs)], label = "Linear", xlabel = "x", ylabel = "Distributed force", legendtitle = "Strain measure", legend = :outerright)
plot(xs, [x[1] for x in lin_force.(xs)],
	label = "Linear",
	xlabel = "Position along the bar",
	ylabel = "Distributed axial force",
	legendtitle = "Strain measure",
	legend = :outerright)
plot!(xs, [x[1] for x in nonlin_force.(xs)], label = "Nonlinear")

savefig(replace(results_file, "results.json" => "forces.tex"))


# strain_limit = 0.02;
strain_limit = maximum(x[1] for x in nonlin_force.(xs)) * linear_problem.length / (3 * bar_E * linear_problem.area);

dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)

linear_problem = fixedBarproblem1D(
	bar_length,
	area,
	lin_force,
	num_ele,
	0.0;
	right_fixed = true,
)

nonlinear_problem = fixedBarproblem1D(
	bar_length,
	area,
	nonlin_force,
	num_ele,
	1.0;
	right_fixed = true,
)
lin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-8)
nonlin_result = Datasolver.greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-8)
plot_results(lin_result, dataset = dataset, title = "Linear solver")
plot_results(nonlin_result, dataset = dataset, title = "Nonlinear solver")
# result = Datasolver.greedyLocalSearchSolverNonLinearBar(linear_problem, dataset)
@show lin_result.cost[end]
@show nonlin_result.cost[end]
# result = NLP_solver(linear_problem, dataset, use_L1_norm = false, worklimit = 100)

xs = [p[1] for p in lin_result.Φ]
# α = f_0 / (bar_E * area * (n * pi / bar_length)^2)
# u_ans = [α * sin(n * pi * x / bar_length) for x in xs]
u_ans = [β * sin(pi * x / bar_length) for x in xs]
@show L2_lin = round(calc_reldiff(xs, Datasolver.get_final(lin_result).u, u_ans), digits = 4)
@show L2_nonlin = round(calc_reldiff(xs, Datasolver.get_final(nonlin_result).u, u_ans), digits = 4)


plot(xs, u_ans, label = "Analytical solution", color = :black, linewidth = 2)
plot!(xs, get_final(lin_result).u, label = "Linear solver", color = :blue, linewidth = 2)
plot!(xs, get_final(nonlin_result).u, label = "NonLinear solver", color = :red, linewidth = 2)

################# Linear #####################
N_datapoints = [2^n for n in 2:10]
N_elements = [2^n for n in 2:10]
analytical_u = []
results = Vector{NamedTuple}()
for N_d in N_datapoints, N_e in N_elements
	@show N_d, N_e
	local dataset = create_dataset(N_d, x -> bar_E * x, -strain_limit, strain_limit)
	local linear_problem = fixedBarproblem1D(
		bar_length,
		area,
		lin_force,
		N_e,
		0.0;
		right_fixed = true,
	)
	local result = greedyLocalSearchSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-9)
	push!(results, (N_datapoints = N_d, N_elements = N_e, result = result))

	xs = [p[1] for p in result.Φ]
	u_ans = [β * sin(pi * x / bar_length) for x in xs]
	push!(analytical_u, u_ans)
end
Datasolver.convergence_analysis(results, analytical_u)


################# Nonlinear #####################

nonlin_analytical_u = []
nonlin_results = Vector{NamedTuple}()
for N_d in N_datapoints, N_e in N_elements
	@show N_d, N_e
	local dataset = create_dataset(N_d, x -> bar_E * x, -strain_limit, strain_limit)
	local nonlinear_problem = fixedBarproblem1D(
		bar_length,
		area,
		nonlin_force,
		N_e,
		1.0;
		right_fixed = true,
	)
	local result = greedyLocalSearchSolverNonLinearBar(nonlinear_problem, dataset, NR_tol = 1e-9)
	push!(nonlin_results, (N_datapoints = N_d, N_elements = N_e, result = result))

	xs = [p[1] for p in result.Φ]
	u_ans = [β * sin(pi * x / bar_length) for x in xs]
	push!(nonlin_analytical_u, u_ans)
end
Datasolver.convergence_analysis(nonlin_results, nonlin_analytical_u)


Datasolver.convergence_analysis_dual(results, analytical_u, nonlin_results, nonlin_analytical_u)

savefig(replace(results_file, "results.json" => "convergence_contour.tex"))
uncomment_pgfplotsset_blocks(dirname(results_file))
