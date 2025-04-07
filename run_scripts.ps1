$juliaScripts = @("examples\NLP_num_datapoints.jl", "examples\optimality_without_noise.jl", "examples\optimality_with_noise.jl", "examples\NLP_objective_strain.jl", "examples\NLP_num_elements.jl", "examples\NLP_strain_measure.jl", "examples\NLP_element_datapoints.jl")

foreach ($script in $juliaScripts) {
    Write-Host "Running $script..."
    try {
        & julia $script
    }
    catch {
        Write-Warning "Error occurred while running $script, continuing..."
    }
}
