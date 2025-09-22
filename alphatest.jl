### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
	#! format: off
	quote
		local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
		local el = $(esc(element))
		global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
		el
	end
	#! format: on
end

# ╔═╡ 4833e9a2-200e-405a-8fa3-081cf21f694e
begin
	using Pkg: Pkg
	Pkg.activate()
	Pkg.develop(path = "./")
	Pkg.instantiate()
	Pkg.add("PlutoUI")
	Pkg.add("StaticArrays")
	Pkg.add("Plots")
	using Datasolver
	using PlutoUI
	using LinearAlgebra, SparseArrays, StaticArrays, Statistics, Plots


end

# ╔═╡ 5ab9ae5b-5d7d-4194-9fad-0931143c30b0
using Printf


# ╔═╡ 90e30260-f295-11ef-1a66-210bdbdb6df1
@bind α Slider(0:0.01:1, default = 1)

# ╔═╡ 95ff341e-47b6-421a-8297-a34b39342a23
md"#### α = $α"

# ╔═╡ e16664bc-2962-4cf0-ad65-92f923d57b07
begin
	NR_steps = 50
	load_steps = 50
end;

# ╔═╡ 362b3733-a4e5-4d46-8887-62756af305c3
md"#### Select dims"

# ╔═╡ 7dfc0397-6aa3-4250-9401-2fb294262a55
@bind dims Select([1, 2, 3])

# ╔═╡ bb367180-19bc-4c1b-9aea-2b61fa99bd86
begin
	bar_area = 0.1     # [m^2] - cross-sectional area of the bar

	bar_len = 9.0      # [m]   - initial length of the bar
	bar_E = 1e5        # [Pa]  - assumed Young_modulus
	num_ele = 10       # [-]   - number of elements
	numDataPts = 201   # [-]   - number of data points
	if dims == 1
		bar_distF = [1.8e2]     # [N]   - constant uniform distributed load
	else
		bar_distF = [1.8e2, 0.0]   # [N]   - constant uniform distributed load
	end

end;

# ╔═╡ af82bef9-80da-4dce-a861-9b31668ffb9d
begin
	strain_limit = 2 * norm(bar_distF) * bar_len / (bar_E * bar_area)
	dataset = dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
	SE = hcat(dataset.E, dataset.S)
	costFunc_ele = (e, s) -> 0.5 * (dataset.C * e^2 + 1 / dataset.C * s^2)

end;

# ╔═╡ 7352fc24-b76d-4303-9f32-f3b52b4c3754
begin

	# node vector
	num_node = num_ele + 1
	# node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
	if dims == 1
		node_vector = [[x] for x in LinRange(0.0, bar_len, num_node)]
		constrained_dofs = get_constrained_dofs([(1, 1)], num_ele, dims)
	else
		node_vector = [[x, 0.0] for x in LinRange(0.0, bar_len, num_node)]
		constrained_dofs = get_constrained_dofs([(1, 1), (1, 2), (num_node, 1), (num_node, 2)], num_ele, dims)
	end
end;

# ╔═╡ 19ab3abe-ac23-4991-b61b-1b27ab937aa3

# # solving
results = directSolverNonLinearBar(
	node_vector,
	constrained_dofs,
	SE,
	costFunc_ele,
	bar_area,
	bar_distF,
	num_ele,
	dataset.C;
	NR_max_iter = NR_steps,
	NR_num_load_step = load_steps,
	alpha = α,
	random_init_data = false,
);



# ╔═╡ 9e5ece3b-abb2-4700-bb32-f49bb7ca97e8

plot_results(results, dataset = dataset)

# ╔═╡ d47816b5-a4d2-419e-812a-4f5ae19910a4
begin
	zero_alpha_results = directSolverNonLinearBar(
		node_vector,
		constrained_dofs,
		SE,
		costFunc_ele,
		bar_area,
		bar_distF,
		num_ele,
		dataset.C;
		NR_max_iter = NR_steps,
		NR_num_load_step = load_steps,
		alpha = 0.0,
		random_init_data = false,
	)

	one_alpha_results = directSolverNonLinearBar(
		node_vector,
		constrained_dofs,
		SE,
		costFunc_ele,
		bar_area,
		bar_distF,
		num_ele,
		dataset.C;
		NR_max_iter = NR_steps,
		NR_num_load_step = load_steps,
		alpha = 1.0,
		random_init_data = false,
	)


end;

# ╔═╡ 91ebddf9-4d2e-45f0-80e1-2a47250f2a2f
begin
	labels = ["α = 0", "α = $α", "α = 1"]
	all_results = [zero_alpha_results, results, one_alpha_results]
end;

# ╔═╡ c52af331-9051-4f2f-be64-d9373099a0bd
function plot_result(results::Vector{SolveResults}, field, labels; title = "")
	# Extract x positions from result.Φ
	x_nodes = [norm(p[1]) for p in results[1].Φ]

	# Calculate midpoints between nodes for lambda and mu
	x_midpoints = [norm((x_nodes[i] + x_nodes[i+1]) / 2) for i in 1:length(x_nodes)-1]

	final_results = [get_final(result) for result in results]
	tick_formatter = x -> @sprintf("%.2g", x)
	# Plot e, s, and u at each node
	if field ∈ [:e, :s, :μ]
		p1 = plot(x_midpoints, final_results[1][field], xlabel = "x", ylabel = String(field), title = String(field), marker = :x, yformatter = tick_formatter, label = labels[1])
		for i in 2:length(results)
			plot!(x_midpoints, final_results[i][field], marker = :x, yformatter = tick_formatter, label = labels[i])
		end

	else
		p1 = plot(x_nodes, final_results[1][field], xlabel = "x", ylabel = String(field), title = String(field), marker = :x, yformatter = tick_formatter, label = labels[1])
		for i in 2:length(results)
			plot!(x_nodes, final_results[i][field], marker = :x, yformatter = tick_formatter, label = labels[i])
		end

	end


	p1

end;


# ╔═╡ 15dbd502-ef04-485c-bf6b-0cbe73edc271
plot_result(all_results, :u, labels)

# ╔═╡ a82d9987-94a8-4345-896e-36395328ec48
plot_result(all_results, :e, labels)

# ╔═╡ dcd669ec-b173-4018-9cb0-a9537b30b471
plot_result(all_results, :s, labels)

# ╔═╡ 7a380121-4d7c-49b3-9dac-1b0b47aa0597
plot_result(all_results, :μ, labels)

# ╔═╡ cc21d645-a38f-4060-b578-3f6bb6212e37
plot_result(all_results, :λ, labels)

# ╔═╡ 3ae841c5-eaac-475d-bd70-247533ad62f0


# ╔═╡ Cell order:
# ╠═4833e9a2-200e-405a-8fa3-081cf21f694e
# ╟─95ff341e-47b6-421a-8297-a34b39342a23
# ╟─90e30260-f295-11ef-1a66-210bdbdb6df1
# ╠═e16664bc-2962-4cf0-ad65-92f923d57b07
# ╠═bb367180-19bc-4c1b-9aea-2b61fa99bd86
# ╠═7352fc24-b76d-4303-9f32-f3b52b4c3754
# ╠═af82bef9-80da-4dce-a861-9b31668ffb9d
# ╟─362b3733-a4e5-4d46-8887-62756af305c3
# ╟─7dfc0397-6aa3-4250-9401-2fb294262a55
# ╠═19ab3abe-ac23-4991-b61b-1b27ab937aa3
# ╟─9e5ece3b-abb2-4700-bb32-f49bb7ca97e8
# ╠═d47816b5-a4d2-419e-812a-4f5ae19910a4
# ╟─91ebddf9-4d2e-45f0-80e1-2a47250f2a2f
# ╟─c52af331-9051-4f2f-be64-d9373099a0bd
# ╠═15dbd502-ef04-485c-bf6b-0cbe73edc271
# ╠═a82d9987-94a8-4345-896e-36395328ec48
# ╠═dcd669ec-b173-4018-9cb0-a9537b30b471
# ╠═7a380121-4d7c-49b3-9dac-1b0b47aa0597
# ╠═cc21d645-a38f-4060-b578-3f6bb6212e37
# ╠═3ae841c5-eaac-475d-bd70-247533ad62f0
# ╟─5ab9ae5b-5d7d-4194-9fad-0931143c30b0
