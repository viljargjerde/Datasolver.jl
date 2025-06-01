include("basic_setup.jl")
using Plots
using StatsPlots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using DelaunayTriangulation


results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
results_list = []

datasets = [create_dataset(numDataPts, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = 0.01) for _ in 1:100]

linear_problem, nonlinear_problem = get_problems(4)
scatter(datasets[1].E, datasets[1].S)
# Define shorthand functions for the solvers
NLP(lin_problem, random_init, data) = NLP_solver(lin_problem ? linear_problem : nonlinear_problem, data; use_L1_norm = false, random_init_data = random_init, worklimit = 200.0, parameter_file = "NLP_params.prm")
directSolver(lin_problem, random_init, data) = directSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init)
greedySearch(lin_problem, random_init, data) = Datasolver.greedyLocalSearchSolverNonLinearBar(lin_problem ? linear_problem : nonlinear_problem, data; random_init_data = random_init, search_iters = 1000)

# solvers = ["ADM", "GO-ADM"]
solvers = ["ADM", "GO-ADM", "MINLP"]
num_exps = length(datasets) * length(solvers) * 2 * 2


if isfile(results_file)
	println("Loading existing results from file...")
	results_list = JSON.parsefile(results_file)
else
	println("Running solver and storing results...")
	results_list = Vector{Dict}(undef, num_exps)
	Threads.@threads for i in eachindex(datasets)
		di = 0
		for lin_problem in [true, false]
			for random_init in [true, false]
				for solver in solvers
					@show solver
					t1 = time()
					if solver == "ADM"
						result = directSolver(lin_problem, random_init, datasets[i])
					elseif solver == "GO-ADM"
						result = greedySearch(lin_problem, random_init, datasets[i])
					elseif solver == "MINLP"
						result = NLP(lin_problem, random_init, datasets[i])
					end

					t2 = time()
					@show result.cost[end]
					results_list[length(solvers)*4*(i-1)+1+di] =
						Dict(
							"Solver" => solver,
							"Initialization" => random_init ? "Random" : "Nullspace",
							"Problem" => lin_problem ? "Linear" : "Nonlinear",
							"Cost" => result.cost[end],
							"Total solve time" => t2 - t1,
							"Model solve time" => result.solvetime[2],
							"Work" => result.solvetime[1],
							"Result" => result,
							"Dataset" => i,
							"S" => datasets[i].S,
							"E" => datasets[i].E,
						)
					di += 1
				end
			end
		end
	end
	# Save results to file
	open(results_file, "w") do f
		JSON.print(f, results_list)
	end
end

df = DataFrame(Dict.(results_list))
combos = unique(select(df, [:Initialization, :Problem]))

function latex_sci(x::Real; digits = 1)
	if x == 0
		return "0"
	end
	exponent = floor(Int, log10(abs(x)))
	base = round(x / 10^exponent, digits = digits)
	return "\$$(base) \\times 10^{$(exponent)}\$"
end
# Create 4 subplots
plot_list = []
row_titles = []
col_titles = []

ymin = minimum(df.Cost)
ymax = maximum(df.Cost)
yticks = range(ymin, ymax; length = 5)


for row in eachrow(combos)
	init = row.Initialization
	prob = row.Problem
	push!(row_titles, "$(prob)")
	push!(col_titles, "$(init)")
	df_sub = filter(r -> r.Initialization == init && r.Problem == prob, df)
	p = @df df_sub boxplot(:Solver, :Cost, legend = false, yticks = yticks, ylim = (ymin, ymax), yformatter = latex_sci)
	println("$prob $init")
	println("ADM: $(sum(filter(r -> r.Solver == "ADM", df_sub).Cost))")
	println("Search: $(sum(filter(r -> r.Solver == "Search", df_sub).Cost))")
	println("GO-ADM: $(sum(filter(r -> r.Solver == "GO-ADM", df_sub).Cost))")
	println("MINLP: $(sum(filter(r -> r.Solver == "MINLP", df_sub).Cost))")
	println()
	push!(plot_list, p)
end


#####################
df_sub = filter(r -> r.Initialization == "Nullspace" && r.Problem == "Nonlinear", df)
df_solver = filter(r -> r.Solver == "ADM", df_sub)
histogram(df_solver.Cost, normalize = :pdf, alpha = 0.5)
@show size(df_sub)
df_solver = filter(r -> r.Solver == "MINLP", df_sub)
histogram!(df_solver.Cost, normalize = :pdf, alpha = 0.5, color = :black)
@show size(df_sub)
df_solver = filter(r -> r.Solver == "GO-ADM", df_sub)
histogram!(df_solver.Cost, normalize = :pdf, alpha = 0.5)
@show size(df_sub)

row_titles = unique(row_titles)
col_titles = unique(col_titles)

plt = plot(plot_list[1:4]..., layout = (2, 2), size = (400, 600))
title!(plt[1], col_titles[1], titlefont = font(12))
title!(plt[2], col_titles[2], titlefont = font(12))

# Add row titles (on y-axis)
ylabel!(plt[1], row_titles[1])
ylabel!(plt[3], row_titles[2])

savefig(plt, replace(results_file, "results.json" => "boxplot.tex"))



############# Table ###############
grouped = groupby(df, [:Initialization, :Problem, :Solver])

stats_df = combine(grouped,
	:Cost => mean => :mean,
	:Cost => median => :median,
	:Cost => std => :std,
	:Cost => minimum => :min,
	:Cost => maximum => :max,
)


table = select(filter(row -> row.Initialization == "Nullspace" && row.Problem == "Linear", stats_df), Not([:Initialization, :Problem]))
pretty_table(table, show_subheader = false)
open(replace(results_file, "results.json" => "nullspace-linear.tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
end

table = select(filter(row -> row.Initialization == "Nullspace" && row.Problem == "Nonlinear", stats_df), Not([:Initialization, :Problem]))
pretty_table(table, show_subheader = false)
open(replace(results_file, "results.json" => "nullspace-nonlinear.tex"), "w") do f
	pretty_table(f, table, backend = Val(:latex), show_subheader = false)
end

#############       ###############

df
df_direct = filter(row -> row.Solver == "ADM" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
df_search = filter(row -> row.Solver == "Search" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
df_greedy = filter(row -> row.Solver == "GO-ADM" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
df_NLP = filter(row -> row.Solver == "MINLP" && row.Problem == "Nonlinear" && row.Initialization == "Nullspace", df)
sort!(df_direct, :Dataset)
sort!(df_search, :Dataset)
sort!(df_greedy, :Dataset)
sort!(df_NLP, :Dataset)
cost_diff = abs.(df_direct.Cost .- df_greedy.Cost)
# cost_diff = abs.(df_direct.Cost .- df_NLP.Cost)
max_diff_idx = argmax(cost_diff)
DE = Vector{Float64}(df_direct[max_diff_idx, :E])
DS = Vector{Float64}(df_direct[max_diff_idx, :S])

### For dataset plot ### 
scatter(DE, DS, label = nothing, xlabel = "Strain [-]", ylabel = "Stress [MPa]", legend = :topleft, marker = :square, line = :dash, markercolor = paired_colors[1])
savefig(replace(results_file, "results.json" => "dataset.tex"))
p = scatter(DE, DS, label = "Data", xlabel = "Strain [-]", ylabel = "Stress [MPa]", legend = :topleft, marker = :square, line = :dash, markercolor = paired_colors[1])
s = Datasolver.get_initialization_s(linear_problem)
best_idxs = Datasolver.find_closest_idx(DS, s)
S = DS[best_idxs]
E = DE[best_idxs]
scatter!(E, S, label = "Initialization", marker = :square, markercolor = paired_colors[5], line = :dash)

r1 = df_direct[max_diff_idx, :Result]
scatter!(r1["e"][end], r1["s"][end], label = "ADM", marker = :diamond, markercolor = paired_colors[9], line = :dash)

r3 = df_greedy[max_diff_idx, :Result]
scatter!(r3["e"][end], r3["s"][end], label = "GO-ADM", marker = :circle, markercolor = paired_colors[12], line = :dash)

r2 = df_NLP[max_diff_idx, :Result]
scatter!(r2["e"][end], r2["s"][end], label = "MINLP", marker = :star, markersize = 2, markercolor = paired_colors[9], line = :dash)




function anisotropic_transform(points, C)
	[(√C * x, y / √C) for (x, y) in points]
end

function inverse_transform(points, C)
	[(x / √C, y * √C) for (x, y) in points]
end

# --- Get anisotropy coefficient and prepare points ---
C = dataset.C
points_orig = collect(zip(DE, DS))
points_aniso = anisotropic_transform(points_orig, C)

# --- Triangulate and compute Voronoi tessellation ---
tri = triangulate(points_aniso)
vorn = voronoi(tri)

# --- Compute bounding box from E, S with padding ---
padding = 0.05  # 5% padding
xmin, xmax = extrema(DE)
ymin, ymax = extrema(DS)
xrange, yrange = xmax - xmin, ymax - ymin
bbox = (
	xmin - padding * xrange,
	xmax + padding * xrange,
	ymin - padding * yrange,
	ymax + padding * yrange,
)
xmin_a, xmax_a = √C * bbox[1], √C * bbox[2]
ymin_a, ymax_a = bbox[3] / √C, bbox[4] / √C
bbox_aniso = (xmin_a, xmax_a, ymin_a, ymax_a)
# --- Plot Voronoi cells ---
for i in each_polygon_index(vorn)
	# coords = i in DelaunayTriangulation.get_unbounded_polygons(vorn) ?
	# 		 get_polygon_coordinates(vorn, i, bbox_aniso) :
	# 		 get_polygon_coordinates(vorn, i)

	coords = get_polygon_coordinates(vorn, i, bbox_aniso)
	coords_orig = inverse_transform(coords, C)
	x, y = first.(coords_orig), last.(coords_orig)
	push!(x, x[1])  # close polygon
	push!(y, y[1])
	plot!(x, y, seriestype = :shape, linealpha = 0.15, fillalpha = 0.00, linecolor = :gray, label = nothing)
end

p

savefig(replace(results_file, "results.json" => "example-solution.tex"))




uncomment_pgfplotsset_blocks(dirname(results_file))



