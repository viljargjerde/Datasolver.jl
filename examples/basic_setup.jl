using LinearAlgebra
using Datasolver

# inputs
bar_length = 1.0      # [m]   - initial length of the bar
area = 0.5     # [m^2] - cross-sectional area of the bar
force = x -> [1.8e2]  # [N]   - constant uniform distributed load
num_ele = 6       # [-]   - number of elements
function get_problems(num_ele)
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
	linear_problem, nonlinear_problem

end

linear_problem, nonlinear_problem = get_problems(num_ele)

bar_E = 1e3;        # [Pa]  - assumed Young_modulus
numDataPts = 21;   # [-]   - number of data points, odd number to ensure zero strain is in the dataset


# Generate data 
dist_F = linear_problem.force(nothing)
strain_limit = norm(dist_F) * linear_problem.length / (bar_E * linear_problem.area);
if norm(linear_problem.force(1)) <= 1e-6
	strain_limit += 1.0
end
dataset = create_dataset(numDataPts, x -> bar_E * x, -strain_limit, strain_limit)


function check_similarity(results1, results2)
	is_same = true
	fields = ["e",
		"E",
		"s",
		"S",
		"u",
		"cost"]
	for field in fields
		if !(results1[field] ≈ results2[field])
			is_same = false
			break
		end
	end
	is_same

end
