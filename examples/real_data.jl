using CSV
using DataFrames
using XLSX
using Plots
using Datasolver
using Revise
filename = "C:/Users/viljar/OneDrive - University of Bergen/Desktop/dyn_stiffness_nylon.xlsx"
gr()

# area = 140.0 # mm^2
area = (140.0 / 2)^2 * π  # mm^2
L_0 = 2076.0
headers = XLSX.readdata(filename, "Sheet1", "B1:I1")[1, :]
data = XLSX.readdata(filename, "Sheet1", "B5:I78960")

df = DataFrame(data, Symbol.(headers))

df = select!(df, ["Time", "Force", "Extensometer [mm]", "% MBL"])


df.Time = df.Time .- df.Time[1]
df.strain = df[!, Symbol("Extensometer [mm]")] ./ L_0

df.stress = df.Force .* 1000 ./ area  # Stress in MPa

# filtered_df = filter(row -> row["% MBL"] < 0.4 && row["stress"] < 50, df)
# filtered_df = filter(row -> row["% MBL"] < 0.4 && row.Time < 177 * 60 , df)
filtered_df = filter(row -> row["% MBL"] < 0.4 && row.Time > 112 * 60 && row.Time < 114 * 60, df)
# filtered_df = filter(row -> row["% MBL"] < 0.4 && row.Time > 2 * 60 && row.Time < 179 * 60, df)



scatter(
	filtered_df.strain,
	filtered_df.stress,
	markersize = 1,
	alpha = 1,
	zcolor = filtered_df.Time .- filtered_df.Time[1],
	label = nothing,
	colorbar = true,  # show colorbar
	xlabel = "Strain [-]",
	ylabel = "Stress [MPa]",
	title = "Stress-Strain colored by Time",
)
# df = DataFrame(XLSX.readtable(sheet)...)

# dataset = Dataset(df.strain, df.stress)
dataset = Dataset(filtered_df.strain, filtered_df.stress)
begin
	nonlinear_problem = fixedBarproblem1D(
		L_0,
		area,
		# x -> [100 * x / L_0],  # [N/mm]   - constant uniform distributed load
		x -> [x > 0.95 * L_0 ? 2500 : 0],  # [N/mm]   - constant uniform distributed load
		# x -> [500e3 / 200 ],  # [N/mm]   - constant uniform distributed load
		16,
		1.0;
		right_fixed = false,
	)
	result_direct = directSolverNonLinearBar(
		nonlinear_problem,
		dataset;
		random_init_data = false,
		verbose = true,
		NR_tol = 1e-9,
	)
	result_greedy = Datasolver.greedyLocalSearchSolverNonLinearBar(
		nonlinear_problem,
		dataset;
		random_init_data = false,
		verbose = true,
		NR_tol = 1e-9,
		search_iters = 1000,
	)
	@show result_direct.cost
	@show result_greedy.cost
	plot_results(result_greedy, dataset = dataset)
end


