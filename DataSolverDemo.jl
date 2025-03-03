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

# ╔═╡ 6b982989-c10d-4b78-af79-df4bc84b6b7a
begin
	import Pkg
	Pkg.activate()
	Pkg.develop(path="./")
	Pkg.instantiate()
	Pkg.add("PlutoUI")
	using PlutoUI
	using Datasolver
	
end

# ╔═╡ 203eea58-f728-4fd9-acc1-c76dd8899473
struct TwoColumn{A, B}
	left::A
	right::B
end

# ╔═╡ 644e61af-6e8a-4ab0-8cb4-0d1e68abd48a
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

# ╔═╡ 0fc2f235-7464-4837-87ab-d834c1622a8d
md"""
# Data driven solver
**Viljar Helgestad Gjerde 04.10.24**

## Contents

* Developed solver package
* Demonstration on new problems
* New solver

"""

# ╔═╡ b785b0eb-cedd-41c8-933d-cd5efa6caa93
md"""
# New package - Datasolver

* Aims to create an interface for easily setting up and solving datadriven problems
* Solver in datasolver.jl
    * Defines types for dataset and results
        * Results stores the history for all relevant variables
    * Defines the datasolve solver
* Several useful functions in utils.jl
    * Setup functions:
        * create_dataset
        * connect\_in\_sequence and create\_$\Phi$_bar
        * Integration functions
    * Visualization functions:
        * plot_dataset with or without results
        * convergence_analysis 
        * plot_results plots everything stored in the results
        * plot_configuration plots nodes with connections
"""

# ╔═╡ 5008014f-e2b9-4b78-81fe-d99c0a5a1ca2
md"""
# How does the new solver work?

1) From the constraint $B^Ts - f = 0 \implies s = B^T \setminus f$
2) Find e from this s
    - We probably don't have this s in the dataset, so we approximate, e.g. linear interpolation
3) From $e - B u = 0$, we find u
"""

# ╔═╡ 2ef5dfad-7057-4cdc-8aa4-00f52b97ccb8
TwoColumn(md"""
## Pros:
* Fast
* Always satisfies constraints exactly 
* Handles noise well
* Works well with few datapoints
* Only needs to solve a single system of equations once
""",md"""
## Cons:
* We no longer select points we have observed in the dataset
* Quality dependent on how the approximation is done 
* No cost measure - hard to evaluate the results
""")

# ╔═╡ 3c306bca-2ed5-479e-b5a0-c660f83f1a46
noise_slider2 = @bind noise2 Slider(0:1e3:1e5,default=0);

# ╔═╡ e18186fb-0beb-403c-b8e1-8f1c48beba74
md"""
# Kanno example

## Setting up the problem
"""

# ╔═╡ 9862ceb4-966c-442f-beb8-d5118a377abc
begin
	L = 3.6 # m
	A = 0.002
	F = 4000 # N
	
	connections = [[1, 3], [3, 5], [3, 4], [5, 6], [2, 4], [4, 6], [1, 4], [2, 3], [3, 6], [4, 5]]
	Φ = [[0, 0], [0, L], [L, 0], [L, L], [2L, 0], [2L, L]]
	
	
	plot_configuration(Φ, connections)
	
end

# ╔═╡ facd1158-a300-4043-8e16-416174acd22c
md"""
## Create the data
"""

# ╔═╡ a46e1a29-9516-4c76-b045-f9bbd444693b
begin
	noise_slider = @bind noise_magnitude Slider(0:1e4:1e6,default=1e5);
	noise_slider
end

# ╔═╡ e627d95e-d26f-4dfa-a085-df67906457f3
md"Noise magnitude = $noise_magnitude"

# ╔═╡ 6219b0fe-0225-4179-b287-f3f862ba2e59
@bind strain_magnitude Slider(0.0001:0.001:0.01,default=0.005)

# ╔═╡ 00406ac5-481f-4240-aa34-9da366fabc71
md"Strain magnitude = $strain_magnitude"

# ╔═╡ 3a8329e8-ef06-46ac-9d99-14aca370614a
begin
	datapoints_slider = @bind N_datapoints Slider(2:500,default=50);
	datapoints_slider
end

# ╔═╡ 9a371c39-1961-4731-838e-5b4699c9189c
md"N datapoints = $N_datapoints"

# ╔═╡ f6e58a57-8460-46b3-984f-210d917a2eab
begin
	dataset = create_dataset(N_datapoints, x -> 5e6 * tanh.(500 .* x), -strain_magnitude, strain_magnitude, noise_magnitude = noise_magnitude)
	plot_dataset(dataset)
	
end

# ╔═╡ fea40ad0-e647-431f-9cbd-bfbc946cebb5
md"""
## Solve the problem and display the choosen points
"""

# ╔═╡ 4113c221-f9d8-47c6-a4a3-8ba60d5e9840
md"""
## This also works using LP with gurobi
"""

# ╔═╡ 36a833a0-8ebe-4358-8b44-73b4a53e2012
let
	fixed_dof = [(1, 1), (1, 2), (2, 1), (2, 2)]
	f = zeros(2 * length(Φ) - 4) 
	f[2] = -F
	f[6] = -F
	_dataset = create_dataset(50, x -> 5e6 * tanh.(500 .* x), -strain_magnitude, strain_magnitude, noise_magnitude = noise_magnitude)
	result = Datasolver.LP_solver(connections, Φ, A, _dataset, f, fixed_dof,verbose=false)
	plot_dataset(dataset, Datasolver.get_final(result))
	
end

# ╔═╡ a3b06c47-cb64-4736-b745-6b76932acbf8
md"""
# Linear bar with continous force
"""

# ╔═╡ 6beabc9c-862c-41f2-ad3d-f7aafcce29c9
md"""
* A bar is discretized into equally sized elements.
* A force at position x is applied according to ``f(x) = f_0 \sin\left(\frac{n \pi x}{L}\right)``.
* Possible to solve analytically:

```math
\begin{align*}
    &-EA u'' = f \\
    &u = \alpha \sin\left(\frac{n \pi x}{L}\right) \\
    &u' = \frac{n \pi}{L} \alpha \cos\left(\frac{n \pi x}{L}\right) \\
    &u'' = \left(\frac{n \pi}{L}\right)^2 \alpha \sin\left(\frac{n \pi x}{L}\right) \\
    &EA \left(\frac{n \pi}{L}\right)^2 \alpha \sin\left(\frac{n \pi x}{L}\right) = f_0 \sin\left(\frac{n \pi x}{L}\right) \\
    &\alpha = \frac{f_0}{EA \left(\frac{n \pi}{L}\right)^2}\\
	& u = \frac{f_0}{EA \left(\frac{n \pi}{L}\right)^2} \sin\left(\frac{n \pi x}{L}\right)
\end{align*}
```
"""


# ╔═╡ 1220d5f1-7134-47df-8892-92c08a0abda0
md"""
## Solution
"""

# ╔═╡ c1e091ac-a26a-410d-9d3e-010747067787
md"N datapoints = $N_datapoints"

# ╔═╡ 6cba6a61-7d12-4147-9947-594c955e8efc
datapoints_slider

# ╔═╡ 27857262-ab5a-4332-b155-e5e0ebda0c20
begin
	elements_slider = @bind N_elements Slider(2:500,default=20);
	elements_slider
end

# ╔═╡ b696a287-ec08-4981-97ae-b6d0ceef2909
md"N elements = $N_elements"

# ╔═╡ 147b3f59-cfaf-4875-8788-4b6fa538eef4
md"""
Noise = $noise2
"""

# ╔═╡ e56eed03-fa28-490a-9cea-1c99ec6dfe0a
noise_slider2

# ╔═╡ 6fcc2fee-2eab-4fef-88dc-854fe7e73fbc
md"""
## Convergence analysis
"""

# ╔═╡ a7e9fd2d-ee2c-4e9c-9976-6f3f5c85833f
md"""
# Different constitutive relations
"""

# ╔═╡ e6a7616a-0135-4239-bd50-ce3d50949926
md"""
## $σ = E_1 ϵ + E_3 ϵ^3$
"""

# ╔═╡ 61a42081-5aaa-497e-ac25-c68f9c331e87
md"N datapoints = $N_datapoints"

# ╔═╡ 5f52fa6c-b870-4faf-9b66-a0ce13d4ba2e
datapoints_slider

# ╔═╡ c6441795-77b2-4b6b-a8e9-2f7987b36b80
md"N elements = $N_elements"

# ╔═╡ 514ba099-4b3a-4b44-af62-29b56578af5c
elements_slider

# ╔═╡ a985ee74-7d1a-43d9-9109-e157d25f7364
md"""
Noise = $noise2
"""

# ╔═╡ 75820a66-759d-452d-a13e-899801d5d428
noise_slider2

# ╔═╡ 2c95f6aa-b992-4cf9-b056-142dc1a73ab9
md"""
## $σ = E \sqrt{|ϵ|} \  sign(ϵ)$
"""

# ╔═╡ 12de5c1e-e1dc-4bf0-b665-b2aafa3190fa
md"N datapoints = $N_datapoints"

# ╔═╡ 98587ba4-46be-4bc9-95fc-891c3d30276a
datapoints_slider

# ╔═╡ 1b22acb9-1ef8-4b96-8c3f-2d4b6466db56
md"N elements = $N_elements"

# ╔═╡ 367e63ec-6025-4c8e-af71-5d4658933efa
elements_slider

# ╔═╡ 7436b325-68c4-4224-bb73-95f0982f1a15
md"""
Noise = $noise2
"""

# ╔═╡ 47f5a159-bc94-422d-8768-503d865bf97f
noise_slider2

# ╔═╡ 8a4613bc-7e5a-43f7-aadf-cef6206ddd89
md"""
##  σ = Choose your own adventure
"""

# ╔═╡ 31ce1fd5-181f-4107-96fe-033edb51e686
σ(ϵ) = 1e6*tanh.(1 .* ϵ)

# ╔═╡ b0d9fcbd-3fdb-4f40-a34c-c9bb844983b9
md"N datapoints = $N_datapoints"

# ╔═╡ 3edb1f2c-e589-4b4b-a8c6-d23642cb507b
datapoints_slider

# ╔═╡ fe72b28c-3547-4955-a516-eaf1fde96116
md"N elements = $N_elements"

# ╔═╡ a4919da0-9895-4f9e-a02d-6a86bbfc6b4a
elements_slider

# ╔═╡ 993913f7-d38a-44b9-8424-7d88c776180c
md"""
Noise = $noise2
"""

# ╔═╡ 31eccca9-6e52-4c61-85f5-9faf72f9ef81
noise_slider2

# ╔═╡ 0c3cdd7e-38b4-4c2a-a772-d8a534b89a18
solver_select = @bind solver Select([Datasolver.datasolve => "Standard solver",Datasolver.my_new_solver=> "New solver",Datasolver.LP_solver => "LP solver"]);

# ╔═╡ d3c95d59-4a51-4cfe-9dc4-13aac6415da1
solver_select

# ╔═╡ c5552c83-3734-4195-a358-50c3779cc5dc
let
	fixed_dof = [(1, 1), (1, 2), (2, 1), (2, 2)]
	f = zeros(2 * length(Φ) - 4) 
	f[2] = -F
	f[6] = -F
	result = solver(connections, Φ, A, dataset, f, fixed_dof,verbose=false)
	plot_dataset(dataset, Datasolver.get_final(result))
	
end

# ╔═╡ 8ac2e487-b99b-426c-9444-32ffb791e43f
solver_select

# ╔═╡ ddae2019-f884-4300-9877-c71c0a5e5d1f
begin
	results_ex2,L2_ex2, dataset_ex2 = let
		L = 1.0
		E = 2e5
		dataset = create_dataset(N_datapoints, x -> E * x, -15.0, 15.0,noise_magnitude=noise2)
		f_0 = 4000.0
		n = 1
		
		
		connections, Φ, f, fixed_dof = setup_1d_bar_problem(N_elements, L, x -> f_0 * sin(n * pi * x / L))
		result = solver(connections, Φ, A, dataset, f, fixed_dof,verbose=false)

		
		xs = [p[1] for p in Φ]
		α = f_0 / (E * A * (n * pi / L)^2)
		u_ans = [α * sin(n * pi * x / L) for x in xs]
		L2 = round(get_rel_diff(xs, [0.0, Datasolver.get_final(result).u..., 0], u_ans),digits=4)
		
		result, L2, dataset
		
	end
	md"Relative L² difference: $(L2_ex2)"
end

# ╔═╡ 9fca0974-1de5-474e-8083-c5f65f90aeeb
plot_results(results_ex2,dataset=dataset_ex2)

# ╔═╡ 3b99edaf-b3dc-4503-95df-f3fd63806fe9
solver_select

# ╔═╡ a964f1e9-324f-46be-9851-140075773524
let
L = 1.0
A = 0.002
E = 2e5
E_3 = 5000.0 #? Supposed to be lower?
dataset = create_dataset(N_datapoints, x -> E * x + E_3 * x^3, -5.0, 5.0;noise_magnitude = noise2)
f_0 = 4000.0
n = 1
connections, Φ, f, fixed_dof = setup_1d_bar_problem(N_elements, L, x -> f_0 * sin(n * pi * x / L))
result = solver(connections, Φ, A, dataset, f, fixed_dof;verbose=false)

# plot_dataset(dataset, Datasolver.get_final(result);title=L"σ = E_1 ϵ + E_3 ϵ^3")
plot_results(result,dataset=dataset)

end

# ╔═╡ 0e910e35-3701-49f9-9b54-1a2c2e066ca8
solver_select

# ╔═╡ 38211d0a-4fe6-4f8a-9034-387b4ec9a05d
let
L = 1.0
A = 0.002
E = 2e5

dataset = create_dataset(N_datapoints, x -> E *sqrt(abs(x))*sign(x), -15.0, 15.0,noise_magnitude=noise2)
f_0 = 4000.0
n = 1
connections, Φ, f, fixed_dof = setup_1d_bar_problem(N_elements, L, x -> f_0 * sin(n * pi * x / L))
result = solver(connections, Φ, A, dataset, f, fixed_dof;verbose=false)


plot_results(result,dataset=dataset)

end

# ╔═╡ 86414094-fd89-45a1-bc23-01d4d07d6d2e
solver_select

# ╔═╡ 1fb87a27-8b68-4b24-9d66-f180a6b8334e
let
L = 1.0
A = 0.002
E = 2e5

dataset = create_dataset(N_datapoints,σ, -15.0, 15.0;noise_magnitude=noise2)
f_0 = 4000.0
n = 1
connections, Φ, f, fixed_dof = setup_1d_bar_problem(N_elements, L, x -> f_0 * sin(n * pi * x / L))
result = solver(connections, Φ, A, dataset, f, fixed_dof,verbose=false)

plot_results(result,dataset=dataset)

end

# ╔═╡ 5bd2eaa0-9930-4207-8f4f-cc7b746c33e2
solver_select_convergence_analysis = @bind solver_convergence Select([Datasolver.datasolve => "Standard solver",Datasolver.my_new_solver=> "New solver",Datasolver.LP_solver => "LP solver"]);

# ╔═╡ d507a121-590f-455e-9925-fade4056bca3
solver_select_convergence_analysis

# ╔═╡ 0a60e7d0-d23f-46c3-860f-d3d22e3a2eb9
let

N_datapoints = [2^n for n in 4:9]
N_elements = [2^n for n in 4:9]
f_0 = 4000.0
n = 1
E = 2e5
analytical_u = []
α = f_0 / (E * A * (n * pi / L)^2)

results = Vector{NamedTuple}()
for N_d in N_datapoints, N_e in N_elements
	local dataset = create_dataset(N_d, x -> E * x, -20.0, 20.0)
	local connections, Φ, f, fixed_dof = setup_1d_bar_problem(N_e, L, x -> f_0 * sin(n * pi * x / L))
	local result = solver_convergence(connections, Φ, A, dataset, f, fixed_dof, verbose = false)
	push!(results, (N_datapoints = N_d, N_elements = N_e, result = result))

	xs = [p[1] for p in Φ]
	u_ans = [α * sin(n * pi * x / L) for x in xs]
	push!(analytical_u, u_ans)
end

convergence_analysis(results, analytical_u)

end

# ╔═╡ Cell order:
# ╠═6b982989-c10d-4b78-af79-df4bc84b6b7a
# ╟─644e61af-6e8a-4ab0-8cb4-0d1e68abd48a
# ╟─203eea58-f728-4fd9-acc1-c76dd8899473
# ╟─0fc2f235-7464-4837-87ab-d834c1622a8d
# ╟─b785b0eb-cedd-41c8-933d-cd5efa6caa93
# ╟─5008014f-e2b9-4b78-81fe-d99c0a5a1ca2
# ╟─2ef5dfad-7057-4cdc-8aa4-00f52b97ccb8
# ╟─3c306bca-2ed5-479e-b5a0-c660f83f1a46
# ╟─e18186fb-0beb-403c-b8e1-8f1c48beba74
# ╠═9862ceb4-966c-442f-beb8-d5118a377abc
# ╟─facd1158-a300-4043-8e16-416174acd22c
# ╟─e627d95e-d26f-4dfa-a085-df67906457f3
# ╟─a46e1a29-9516-4c76-b045-f9bbd444693b
# ╟─00406ac5-481f-4240-aa34-9da366fabc71
# ╟─6219b0fe-0225-4179-b287-f3f862ba2e59
# ╟─9a371c39-1961-4731-838e-5b4699c9189c
# ╠═3a8329e8-ef06-46ac-9d99-14aca370614a
# ╠═f6e58a57-8460-46b3-984f-210d917a2eab
# ╟─fea40ad0-e647-431f-9cbd-bfbc946cebb5
# ╟─d3c95d59-4a51-4cfe-9dc4-13aac6415da1
# ╠═c5552c83-3734-4195-a358-50c3779cc5dc
# ╟─4113c221-f9d8-47c6-a4a3-8ba60d5e9840
# ╠═36a833a0-8ebe-4358-8b44-73b4a53e2012
# ╟─a3b06c47-cb64-4736-b745-6b76932acbf8
# ╟─6beabc9c-862c-41f2-ad3d-f7aafcce29c9
# ╟─1220d5f1-7134-47df-8892-92c08a0abda0
# ╟─c1e091ac-a26a-410d-9d3e-010747067787
# ╟─6cba6a61-7d12-4147-9947-594c955e8efc
# ╟─b696a287-ec08-4981-97ae-b6d0ceef2909
# ╟─27857262-ab5a-4332-b155-e5e0ebda0c20
# ╟─147b3f59-cfaf-4875-8788-4b6fa538eef4
# ╟─e56eed03-fa28-490a-9cea-1c99ec6dfe0a
# ╟─8ac2e487-b99b-426c-9444-32ffb791e43f
# ╠═ddae2019-f884-4300-9877-c71c0a5e5d1f
# ╠═9fca0974-1de5-474e-8083-c5f65f90aeeb
# ╟─6fcc2fee-2eab-4fef-88dc-854fe7e73fbc
# ╟─d507a121-590f-455e-9925-fade4056bca3
# ╠═0a60e7d0-d23f-46c3-860f-d3d22e3a2eb9
# ╟─a7e9fd2d-ee2c-4e9c-9976-6f3f5c85833f
# ╟─e6a7616a-0135-4239-bd50-ce3d50949926
# ╟─61a42081-5aaa-497e-ac25-c68f9c331e87
# ╟─5f52fa6c-b870-4faf-9b66-a0ce13d4ba2e
# ╟─c6441795-77b2-4b6b-a8e9-2f7987b36b80
# ╟─514ba099-4b3a-4b44-af62-29b56578af5c
# ╟─a985ee74-7d1a-43d9-9109-e157d25f7364
# ╟─75820a66-759d-452d-a13e-899801d5d428
# ╟─3b99edaf-b3dc-4503-95df-f3fd63806fe9
# ╟─a964f1e9-324f-46be-9851-140075773524
# ╟─2c95f6aa-b992-4cf9-b056-142dc1a73ab9
# ╟─12de5c1e-e1dc-4bf0-b665-b2aafa3190fa
# ╟─98587ba4-46be-4bc9-95fc-891c3d30276a
# ╟─1b22acb9-1ef8-4b96-8c3f-2d4b6466db56
# ╟─367e63ec-6025-4c8e-af71-5d4658933efa
# ╟─7436b325-68c4-4224-bb73-95f0982f1a15
# ╟─47f5a159-bc94-422d-8768-503d865bf97f
# ╟─0e910e35-3701-49f9-9b54-1a2c2e066ca8
# ╟─38211d0a-4fe6-4f8a-9034-387b4ec9a05d
# ╟─8a4613bc-7e5a-43f7-aadf-cef6206ddd89
# ╠═31ce1fd5-181f-4107-96fe-033edb51e686
# ╟─b0d9fcbd-3fdb-4f40-a34c-c9bb844983b9
# ╠═3edb1f2c-e589-4b4b-a8c6-d23642cb507b
# ╟─fe72b28c-3547-4955-a516-eaf1fde96116
# ╟─a4919da0-9895-4f9e-a02d-6a86bbfc6b4a
# ╟─993913f7-d38a-44b9-8424-7d88c776180c
# ╟─31eccca9-6e52-4c61-85f5-9faf72f9ef81
# ╟─86414094-fd89-45a1-bc23-01d4d07d6d2e
# ╠═1fb87a27-8b68-4b24-9d66-f180a6b8334e
# ╟─0c3cdd7e-38b4-4c2a-a772-d8a534b89a18
# ╟─5bd2eaa0-9930-4207-8f4f-cc7b746c33e2
