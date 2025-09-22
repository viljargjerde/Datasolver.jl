module Datasolver


include("solver.jl")
export Dataset, SolveResults, datasolve, my_new_solver, get_final
include("utils.jl")
export integrate,create_dataset,connect_in_sequence,create_Φ_bar,get_integration_interval,plot_configuration,discretice_1d_force,get_integration_interval,setup_1d_bar_problem,plot_dataset,plot_dataset,get_rel_diff,convergence_analysis,plot_results
end
