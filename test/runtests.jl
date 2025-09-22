using LinearAlgebra
using Revise
using Datasolver
using Test

using Test

tests = [
	"assembly.jl",
	"base.jl",
	"NewtonRaphsonStep.jl",
	"test_utils.jl",
	"test_dataset.jl",
	"test_solver.jl",
	# "non_lin_strain_tests.jl",
]

@testset "Datasolver" begin
	for t in tests
		fp = joinpath(dirname(@__FILE__), t)
		println("$fp ...")
		include(fp)
	end
end; # @testset
