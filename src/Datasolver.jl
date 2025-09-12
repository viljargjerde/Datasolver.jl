module Datasolver

include("dataset.jl")
include("dataproblem.jl")
export Barproblem1D, fixedBarproblem1D, Dataproblem
include("utils.jl")
export Dataset, SolveResults, get_final
include("assembly.jl")
include("solver.jl")
include("LP_solver.jl")
export NLP_solver
export
    create_dataset, plot_dataset, calc_reldiff, convergence_analysis, plot_results
export get_constrained_dofs
export directSolverNonLinearBar, greedyLocalSearchSolverNonLinearBar

end
