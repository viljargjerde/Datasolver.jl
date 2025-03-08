module Datasolver

include("dataset.jl")
export Barproblem, Barproblem1D, fixedBarproblem1D
include("barproblem.jl")
include("utils.jl")
export Dataset, SolveResults, get_final
include("assembly.jl")
include("solver.jl")
include("LP_solver.jl")
export NLP_solver
export
	create_dataset, plot_dataset, get_rel_diff, convergence_analysis, plot_results
export get_constrained_dofs
export directSolverNonLinearBar

end
