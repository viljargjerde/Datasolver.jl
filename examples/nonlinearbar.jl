using Revise
using LinearAlgebra, SparseArrays, StaticArrays, Plots
using Datasolver

#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law

# inputs
length = 1.0      # [m]   - initial length of the bar
area = 0.5     # [m^2] - cross-sectional area of the bar
force = x -> [1.8e2]  # [N]   - constant uniform distributed load
num_ele = 8       # [-]   - number of elements
alpha = 1.0

problem = fixedBarproblem1D(
	length,
	area,
	force,
	num_ele,
	alpha;
	right_fixed = true,
)
bar_E = 1e3;        # [Pa]  - assumed Young_modulus
numDataPts = 201;   # [-]   - number of data points, odd number to ensure zero strain is in the dataset


# Generate data 
dist_F = problem.force(nothing)
strain_limit = 2 * norm(dist_F) * problem.length / (bar_E * problem.area);
if norm(problem.force(1)) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E / 4 * tanh.(10x), -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)
# dataset = create_dataset(numDataPts, x -> bar_E * tanh.(10x), -strain_limit, strain_limit)

# node vector
num_node = num_ele + 1;
node_vector = [[x] for x in LinRange(0.0, problem.length, num_node)]
constrained_dofs = get_constrained_dofs([(1, 1), (num_node, 1)], problem.num_ele, problem.dims)

# # solving
results = directSolverNonLinearBar(
	problem,
	dataset;
	NR_max_iter = 100,
	NR_num_load_step = 100,
	random_init_data = false,
	DD_max_iter = 10,
	NR_damping = 0.5,
);

plot(results.u)
plot_results(results, dataset = dataset)

