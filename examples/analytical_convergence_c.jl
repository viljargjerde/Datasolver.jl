include("basic_setup.jl")
using Datasolver
using JSON
using Statistics
using Plots
using LaTeXStrings
results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
mkpath(dirname(results_file))
# Problem parameters
β = 0.01
bar_length = 1.0
area = 1.0
bar_E = 1.0e3
numDataPts = 32
num_ele = 4

force_func(α, x) =
	[bar_E * area * ((((1 / 2 * β) * pi^(2)) * sin((pi * x / bar_length))) * (((3 * α^(2)) * (((β * pi) * cos((pi * x / bar_length)) / bar_length))^(2)) + ((((6 * α) * β) * pi) * cos((pi * x / bar_length)) / bar_length) + 2) / bar_length^(2))]
lin_force = x -> force_func(0.0, x)
nonlin_force = x -> force_func(1.0, x)

xs = 0:0.01:bar_length
strain_limit = maximum(x[1] for x in lin_force.(xs)) * bar_length / (3 * bar_E * area)
dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit, noise_magnitude = 0.01)

# Convergence parameter ranges
element_range = [2^n for n in 2:6]
datapoint_range = [2^n for n in 2:8]
c_factors = [0.9, 1.0, 1.1]


# Load previous results
if isfile(results_file)
	println("Loading existing results...")
	results_list = JSON.parsefile(results_file)
else
	results_list = []
end

println("Running solver and storing results...")
for num_data_pts in datapoint_range
	dataset = create_dataset(num_data_pts, x -> bar_E * x, -strain_limit, strain_limit, noise_magnitude = 0.01)
	base_c = dataset.C
	for c_fac in c_factors
		dataset = Datasolver.Dataset(dataset.E, dataset.S, c_fac * base_c)

		for num_ele in element_range
			# Skip if result already stored
			if any(r -> r["Datapoints"] == num_data_pts && r["Elements"] == num_ele && r["c"] ≈ c_fac, results_list)
				continue
			end

			println("Running: c=$(round(c_fac, sigdigits=3)), data=$num_data_pts, ele=$num_ele")
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
				0.0;
				right_fixed = true,
			)

			t1 = time()
			result = greedyLocalSearchSolverNonLinearBar(linear_problem, dataset, NR_tol = 1e-9)
			t2 = time()

			xs = [p[1] for p in result.Φ]
			u_ref = [β * sin(pi * x / bar_length) for x in xs]
			u_num = get_final(result).u
			L2_error = calc_reldiff(xs, u_num, u_ref)

			push!(results_list, Dict(
				"Datapoints" => num_data_pts,
				"Elements" => num_ele,
				"c" => c_fac,
				"Mean Time" => t2 - t1,
				"L2 Error" => L2_error,
				"Result" => result,
			))

			# Save after each run
			open(results_file, "w") do f
				JSON.print(f, results_list)
			end
		end
	end
end

# === Visualization ===
using DataFrames
df = DataFrame(results_list)
grouped = groupby(df, [:c])



plots = []

function get_heatmap_data(g)
	pivoted = unstack(g, :Elements, :Datapoints, "L2 Error")
	heatmap_data = Matrix(select(pivoted, Not(:Elements)))
	return heatmap_data
end

heatmap_data = [get_heatmap_data(g) for g in grouped]
base = heatmap_data[2]
vmin = min(minimum((heatmap_data[1] .- base) ./ base), minimum((heatmap_data[3] .- base) ./ base))
vmax = max(maximum((heatmap_data[1] .- base) ./ base), maximum((heatmap_data[3] .- base) ./ base))

for i in [1, 3]
	c_val = unique(grouped[i].c)[1]
	@show maximum((heatmap_data[i] .- base) ./ base)
	@show mean((heatmap_data[i] .- base) ./ base)
	show_xticks = (i == length(grouped))  # only bottom plot shows x-axis
	p = heatmap(
		datapoint_range,
		element_range,
		(heatmap_data[i] .- base) ./ base, ;
		xlabel = show_xticks ? "Data points" : "",
		ylabel = "Elements",
		colorbar_title = L"$\Delta L^2$ error",
		title = L"%$c_val $\times$ c",
		c = :viridis,
		clims = (vmin, vmax),
		xscale = :log2,
		yscale = :log2,
		xticks = datapoint_range,  # hide xticks if not last
		yticks = element_range,
		grid = false,
	)
	push!(plots, p)
end



# Combine all plots into one figure with consistent color scale
combined_plot = plot(plots...; layout = (2, 1), size = (400, 600), link = :x)


for i in 1:size(grouped[1], 1)
	if !(grouped[1][i, "L2 Error"] ≈ grouped[2][i, "L2 Error"]) || !(grouped[1][i, "L2 Error"] ≈ grouped[3][i, "L2 Error"])
		@show grouped[1][i, "Elements"], grouped[1][i, "Datapoints"]
	end
end

uncomment_pgfplotsset_blocks(dirname(results_file))
