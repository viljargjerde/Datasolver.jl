### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 2d1cbd50-002a-11f0-3513-0fbcdf4e4cb3
# ╠═╡ show_logs = false
begin
	using Pkg: Pkg
	Pkg.activate()
	Pkg.develop(path = "./")
	Pkg.instantiate()
	Pkg.add("PlutoUI")
	Pkg.add("PrettyTables")
	Pkg.add("DataFrames")
	Pkg.add("JSON")
	using PlutoUI
	using Datasolver
	using PrettyTables
	using DataFrames
	using JSON
	using Statistics


end

# ╔═╡ 4a2f447c-730a-4553-8535-346eb4ea5718
struct TwoColumn{A, B}
	left::A
	right::B
end

# ╔═╡ f7a0dac5-319e-431b-943b-9dc1ae93df14
function Base.show(io, mime::MIME"text/html", tc::TwoColumn)
	write(io,
		"""
		<div style="display: flex;">
			<div style="flex: 50%;">
		""")
	show(io, mime, tc.left)
	write(io,
		"""
			</div>
			<div style="flex: 50%;">
		""")
	show(io, mime, tc.right)
	write(io,
		"""
			</div>
		</div>
	""")
end

# ╔═╡ f9563f33-941d-44d7-a65c-c4a63d319f74
function load_table(file_path, name_symbol, name_string)

	results_list = JSON.parsefile(file_path)

	# Convert results to DataFrame
	df = DataFrame(Dict.(results_list))

	# Compute statistics
	grouped = combine(groupby(df, [name_symbol, :random_init]),
		:solve_time => median => :median_solve_time,
		:solve_time => (x -> quantile(x, 0.25)) => :q25,
		:solve_time => (x -> quantile(x, 0.75)) => :q75,
	)

	table = unstack(grouped, name_symbol, :random_init, :median_solve_time)
	rename!(table, Dict(name_symbol => name_string))
	rename!(table, Dict("true" => "Random initialization"))
	rename!(table, Dict("false" => "Non-Random initialization"))

	select!(table, name_string, "Random initialization", "Non-Random initialization")
	table

end;

# ╔═╡ fa7b21b6-defa-4be5-86a9-31739d94271f
function render_table(table)
	pretty_table(HTML, table, show_subheader = false)
end;

# ╔═╡ 68b52a1c-c737-427a-9dd6-8ae7b84ee0ef
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
		padding-left: max(160px, 20%);
		padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 00950cb7-0f16-4da5-85a6-3f667c4d7d53
md"""
# Contents

* NLP
	- Problem setting
	- Solve time vs number of datapoints
	- Solve time vs number of elements 
	- Cost of Nonlinear strain measure  
	- L1 objective function vs L2 
* Alternating direction solver
	- TODO
"""

# ╔═╡ 76184d6d-ae1b-4625-a840-41f12bbef6f4
md"# NLP"

# ╔═╡ 9c607ac5-56ba-4028-9343-41c893d2f8f1
md"""
# Problem setting (unless otherwise specified)

-  $L = 1$
-  $A = 0.5$
-  $E = 10^3$
-  Constant distributed force $F = 1.8 \times 10^2$
-  Number of elements = 6
-  Number of data points = 21
-  1D bar with both ends fixed
-  Linear strain measure
-  All timing results are given as the median with 25th and 75th percentile errorbars based on 10 repetitions
"""


# ╔═╡ 4087dee9-9ea7-479a-be5c-1468a882084c
md"""
# Scaling wrt. number of datapoints
* Initializing the binary variables for choosen data points gives ~ 2x-4x speedup 
* Initializing the continous variables hurts performance
* Roughly linear scaling, yet not consistent
"""

# ╔═╡ efcc6513-2f43-4a25-956d-5177c24f48b7
# ╠═╡ show_logs = false
TwoColumn(md"""
$(LocalResource("examples/NLP_solver_comparisons/initialization_num_datapoints/results.png"))""", render_table(load_table("examples/NLP_solver_comparisons/initialization_num_datapoints/results.JSON", :num_data_pts, "Number of datapoints")))

# ╔═╡ de043c5b-34cb-4f02-89e0-dd74ba3a3eba
md"""
# Solve time vs number of elements
* Solve time quickly explodes
* Initializing helps somewhat
"""

# ╔═╡ ff3ebf2f-3aa3-41c6-a641-6c8e1228390a
# ╠═╡ show_logs = false
TwoColumn(md"""
$(LocalResource("examples/NLP_solver_comparisons/initialization_num_elements/results.png"))""", render_table(load_table("examples/NLP_solver_comparisons/initialization_num_elements/results.JSON", :num_ele, "Number of elements")))

# ╔═╡ 232534e3-2a18-434a-900a-f3e5ab1bf6e6
md"""
# Cost of Nonlinear strain measure
* Number of elements reduced to 5
* Feasiblity of solving reduced drastically
* 10 elements impossible
"""

# ╔═╡ 61a1ea0f-d3bc-4c8b-95c9-0d5659e728b5
let
	table = load_table("examples/NLP_solver_comparisons/initialization_strain_measure/results.JSON", :is_non_linear, "Strain measure")
	table."Strain measure" = replace(table."Strain measure", false => "linear", true => "Nonlinear")
	render_table(table)
end

# ╔═╡ 6899100c-ed35-4c1b-89b9-9d1731dbba7a
md"""
# L1 objective function vs L2
* Number of elements increased to 20
* L1 faster than L2
* Choosen data points almost always the same
* strain and stress differ slightly
"""

# ╔═╡ ab450559-5e3c-439f-b8b3-3d6407e6bbc7
let
	table = load_table("examples/NLP_solver_comparisons/initialization_objective_function/results.JSON", :is_L1, "Objective function")
	table."Objective function" = replace(table."Objective function", true => "L2", false => "L1")
	render_table(table)
end

# ╔═╡ Cell order:
# ╟─2d1cbd50-002a-11f0-3513-0fbcdf4e4cb3
# ╟─f7a0dac5-319e-431b-943b-9dc1ae93df14
# ╟─4a2f447c-730a-4553-8535-346eb4ea5718
# ╟─f9563f33-941d-44d7-a65c-c4a63d319f74
# ╟─fa7b21b6-defa-4be5-86a9-31739d94271f
# ╟─68b52a1c-c737-427a-9dd6-8ae7b84ee0ef
# ╟─00950cb7-0f16-4da5-85a6-3f667c4d7d53
# ╟─76184d6d-ae1b-4625-a840-41f12bbef6f4
# ╟─9c607ac5-56ba-4028-9343-41c893d2f8f1
# ╟─4087dee9-9ea7-479a-be5c-1468a882084c
# ╟─efcc6513-2f43-4a25-956d-5177c24f48b7
# ╟─de043c5b-34cb-4f02-89e0-dd74ba3a3eba
# ╟─ff3ebf2f-3aa3-41c6-a641-6c8e1228390a
# ╟─232534e3-2a18-434a-900a-f3e5ab1bf6e6
# ╟─61a1ea0f-d3bc-4c8b-95c9-0d5659e728b5
# ╟─6899100c-ed35-4c1b-89b9-9d1731dbba7a
# ╟─ab450559-5e3c-439f-b8b3-3d6407e6bbc7
