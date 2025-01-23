module Datasolver


include("solver.jl")
export Dataset, SolveResults, datasolve, my_new_solver, get_final, AlgebraicFunction
include("utils.jl")
export integrate,
	create_dataset, connect_in_sequence, create_Φ_bar, get_integration_interval, plot_configuration, discretice_1d_force, get_integration_interval, setup_1d_bar_problem, plot_dataset, plot_dataset, get_rel_diff, convergence_analysis, plot_results
include("DataDrivenNonlinearBar.jl")
export assembleLinearSystemMatrix, assembleRhsLinearBar
export assembleBalanceResidual, assembleLinearizedSystemMatrix
export linearLagrangePolynomials, compute1stDeriv4linearLagrangePolynomials, constantFunctions
export constructBasisFunctionMatrixLinearLagrange, constructBasisFunctionMatrixConstantFuncs
export GaussLegendreQuadRule
export NewtonRaphsonStep
export directSolverLinearBar, directSolverNonLinearBar

end
