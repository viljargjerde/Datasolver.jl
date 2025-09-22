

# Test for remove_dofs(B::Matrix, node_coordinate_pairs) -> Matrix
@testset "remove_dofs Matrix Tests" begin
	B = [1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0; 9.0 10.0 11.0 12.0]
	node_coordinate_pairs = [(2, 1)]

	new_B = Datasolver.remove_dofs(B, node_coordinate_pairs)
	@test new_B == [1.0 2.0 4.0; 5.0 6.0 8.0; 9.0 10.0 12.0]
end



# Test for remove_dofs(f::Vector, node_coordinate_pairs) -> Vector
@testset "remove_dofs Vector Tests" begin
	f = [10.0, 20.0, 30.0, 40.0]
	node_coordinate_pairs = [(2, 1)]
	new_f = Datasolver.remove_dofs(f, node_coordinate_pairs)
	@test new_f == [10.0, 20.0, 40.0]

	node_coordinate_pairs = [(2, 2)]
	new_f = Datasolver.remove_dofs(f, node_coordinate_pairs)
	@test new_f == [10.0, 20.0, 30.0]

	node_coordinate_pairs = [(1, 1), (1, 2), (2, 2)]
	new_f = Datasolver.remove_dofs(f, node_coordinate_pairs)
	@test new_f == [30.0]
end
