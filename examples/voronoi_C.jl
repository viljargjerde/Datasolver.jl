include("basic_setup.jl")
using Plots
using StatsPlots
using DataFrames
using Statistics
using CSV
using JSON
using PrettyTables
using DelaunayTriangulation

# gr() #!
results_file = joinpath("../master_thesis/figures/", splitext(basename(@__FILE__))[1], "results.json")
mkpath(dirname(results_file))


function create_plot(num_data, noise_mag)

	dataset = create_dataset(num_data, x -> bar_E * x, 0.0, 2 * strain_limit, noise_magnitude = noise_mag)
	# datasets = [create_dataset(numDataPts, x -> bar_E * tanh.(20x), -strain_limit, strain_limit, noise_magnitude = 0.2) for _ in 1:10]

	DE, DS = dataset.E, dataset.S
	p = scatter(DE, DS, label = nothing, markersize = 2, xlabel = "Strain [-]", ylabel = "Stress [MPa]")

	function anisotropic_transform(points, C)
		[(√C * x, y / √C) for (x, y) in points]
	end

	function inverse_transform(points, C)
		[(x / √C, y * √C) for (x, y) in points]
	end

	# --- Get anisotropy coefficient and prepare points ---

	for (color, factor) in zip([single_colors[1], single_colors[3], single_colors[5]], [0.9, 1.0, 1.1])
		C = dataset.C * factor
		points_orig = collect(zip(DE, DS))
		points_aniso = anisotropic_transform(points_orig, C)

		# --- Triangulate and compute Voronoi tessellation ---
		@show "Computing Voronoi tessellation for factor $factor"
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
		@show "Drawing Voronoi cells for factor $factor"
		# --- Plot Voronoi cells ---
		# plot!(seriestype = :line, linealpha = 0.25, fillalpha = 0.00, linecolor = color, label = "$factor", legend = :bottomright)
		for i in each_polygon_index(vorn)
			# coords = i in DelaunayTriangulation.get_unbounded_polygons(vorn) ?
			# 		 get_polygon_coordinates(vorn, i, bbox_aniso) :
			# 		 get_polygon_coordinates(vorn, i)

			coords = get_polygon_coordinates(vorn, i, bbox_aniso)
			coords_orig = inverse_transform(coords, C)
			x, y = first.(coords_orig), last.(coords_orig)
			push!(x, x[1])  # close polygon
			push!(y, y[1])
			plot!(x, y, seriestype = :shape, linealpha = 0.35, label = i == 1 ? "$factor × c" : nothing, fillalpha = 0.00, linecolor = color)
		end
	end
	return p
end
p = create_plot(50, 0.05)
savefig(replace(results_file, "results.json" => "50_005.tex"))
p = create_plot(10, 0.01);
savefig(replace(results_file, "results.json" => "10_001.tex"))
p = create_plot(20, 0.05);
p
# p


uncomment_pgfplotsset_blocks(dirname(results_file))
