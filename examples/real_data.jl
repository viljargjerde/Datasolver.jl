using CSV
using DataFrames
using XLSX
using Plots
using DelaunayTriangulation
using Datasolver
using Revise
filename = "C:/Users/viljar/OneDrive - University of Bergen/Desktop/dyn_stiffness_nylon.xlsx"
include("basic_setup.jl")
# gr()
results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
mkpath(dirname(results_file))
dataset = create_dataset(10, x -> 100 * x, 0.0, 2.0, noise_magnitude = 0.01)

function voronoi_plot(dataset)

	DE, DS = dataset.E, dataset.S
	p = scatter(DE, DS, label = nothing, markersize = 2)

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

		coords = get_polygon_coordinates(vorn, i, bbox_aniso)
		coords_orig = inverse_transform(coords, C)
		x, y = first.(coords_orig), last.(coords_orig)
		push!(x, x[1])  # close polygon
		push!(y, y[1])
		plot!(x, y, seriestype = :shape, linealpha = 0.35, fillalpha = 0.00, linecolor = single_colors[2], label = nothing)
	end
	return p
end


# area = 140.0 # mm^2
area = (140.0 / 2)^2 * π  # mm^2
L_0 = 2076.0
headers = XLSX.readdata(filename, "Sheet1", "B1:I1")[1, :]
data = XLSX.readdata(filename, "Sheet1", "B5:I78960")

df = DataFrame(data, Symbol.(headers))

df = select!(df, ["Time", "Force", "Extensometer [mm]", "% MBL"])


df.strain = df[!, Symbol("Extensometer [mm]")] ./ L_0
df.strain = df.strain .- df.strain[argmin(df.strain)]
df.strain = vcat(fill(NaN, 2), df.strain[1:end-2])
df = df[3:end, :]  # Remove first two rows with NaN values
df.stress = df.Force .* 1000 ./ area  # Stress in MPa

df.Time = df.Time .- df.Time[1]

filtered_df = filter(row -> row["% MBL"] < 0.2 && row.Time > 108 * 60 && row.Time < 110 * 60, df)

scatter(
	filtered_df.strain,
	filtered_df.stress,
	markersize = 1.5,
	# alpha = 0.7,
	markerstrokewidth = 0.5,
	zcolor = filtered_df.Time .- filtered_df.Time[1],
	label = nothing,
	colorbar = true,  # show colorbar
	# c = :viridis,  # show colorbar
	colorbar_title = "Time [s]",  # show colorbar
	xlabel = "Strain [-]",
	ylabel = "Stress [MPa]")
savefig(replace(results_file, "results.json" => "several-cycles.tex"))


filtered_df = filter(row -> row["% MBL"] < 0.2 && row.Time > 108 * 60 && row.Time < 108.5 * 60, df)
force_diff = [0.0; diff(filtered_df.Force)]
filtered_df = filtered_df[force_diff.>1, :]
dataset = Dataset(filtered_df.strain, filtered_df.stress)

voronoi_plot(dataset)

begin
	num_ele = 8
	nonlinear_problem = fixedBarproblem1D(
		L_0,
		area,
		x -> [x >= 0.95 * L_0 ? 1200.0 : 30.0],  # [N/mm]   - constant uniform distributed load
		num_ele,
		1.0;
		right_fixed = false,
	)
	result_direct = directSolverNonLinearBar(
		nonlinear_problem,
		dataset;
		random_init_data = false,
		verbose = false,
		NR_tol = 1e-9,
	)
	result_greedy = Datasolver.greedyLocalSearchSolverNonLinearBar(
		nonlinear_problem,
		dataset;
		random_init_data = false,
		verbose = false,
		NR_tol = 1e-9,
		search_iters = 1000,
	)

	@show result_direct.cost
	@show result_greedy.cost
	plot_results(result_greedy, dataset = dataset)
end
greedy_final = get_final(result_greedy)
direct_final = get_final(result_direct)
scatter(
	filtered_df.strain,
	filtered_df.stress,
	markersize = 3,
	label = "Dataset",
	marker = :square,
	xlabel = "Strain [-]",
	ylabel = "Stress [MPa]")
scatter!(greedy_final.e, greedy_final.s, markersize = 3, color = paired_colors[12], label = "GO-ADM solution")
savefig(replace(results_file, "results.json" => "resultscatter.tex"))
xs = collect(0:5:nonlinear_problem.length)
plot(xs ./ nonlinear_problem.length, ([nonlinear_problem.force(x)[1] for x in xs]))

uncomment_pgfplotsset_blocks(dirname(results_file))
