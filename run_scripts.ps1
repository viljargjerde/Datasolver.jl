$juliaScripts = @("examples/analytical_convergence_c.jl", "examples/direct_element_datapoints_bigger.jl", "examples/direct_element_datapoints.jl", "examples/direct_num_datapoints.jl", "examples/direct_num_elements.jl", "examples/direct_strain_measure_iters.jl", "examples/direct_strain_measure.jl", "examples/NLP_element_datapoints.jl", "examples/NLP_num_datapoints.jl", "examples/NLP_num_elements_fixed.jl", "examples/NLP_num_elements.jl", "examples/NLP_objective_function.jl", "examples/NLP_objective_strain.jl", "examples/NLP_strain_measure.jl", "examples/optimality_with_noise_bigger.jl", "examples/optimality_with_noise_heatmap.jl", "examples/optimality_with_noise_tanh.jl", "examples/optimality_with_noise.jl", "examples/optimality_without_noise.jl", "examples/real_data.jl", "examples/S_NLP_improvement_time.jl", "examples/SA_ADM_improvement_time.jl", "examples/SA_ADM_justification_example.jl", "examples/voronoi_C.jl" )

foreach ($script in $juliaScripts) {
    Write-Host "Running $script..."
    try {
        & julia $script
    }
    catch {
        Write-Warning "Error occurred while running $script, continuing..."
    }
}
