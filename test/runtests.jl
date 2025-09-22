using LinearAlgebra
using Revise
using Datasolver
using Test

using Test

tests = [
	"test_assembly.jl",
	"test_base.jl",
	"test_NewtonRaphsonStep.jl",
	"test_utils.jl",
	"test_dataset.jl",
	"test_solver.jl",
]

@testset "Datasolver" begin
	for t in tests
		fp = joinpath(dirname(@__FILE__), t)
		println("$fp ...")
		include(fp)
	end
end; # @testset
