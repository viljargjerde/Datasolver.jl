

struct Dataproblem{T<:Union{Function,AbstractVector}}
    length::Float64
    area::Float64
    force::T
    num_ele::Int64
    num_node::Int64
    alpha::Float64
    connections::Vector{Tuple{Int64,Int64}}
    constrained_dofs::Vector{Int64}
    node_vector::Vector{Vector{Float64}}
    num_quad_pts::Int64
    dims::Int64
end



function Barproblem1D(
    length::Float64,
    area::Float64,
    force::Function,
    num_ele::Int64,
    alpha::Float64,
    constrained_dofs::Vector{Int64};
    node_vector::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    num_quad_pts::Int64=2,
)
    num_node = num_ele + 1
    if isempty(node_vector)
        node_vector = [[x] for x in LinRange(0.0, length, num_node)]
    end
    connections = [(i, i + 1) for i in 1:num_ele]
    return Dataproblem(length, area, force, num_ele, num_node, alpha, connections, constrained_dofs, node_vector, num_quad_pts, 1)
end


function TrussProblem(
    area::Float64,
    force::AbstractVector,
    connections::Vector{Tuple{Int64,Int64}},
    alpha::Float64,
    constrained_dofs::Vector{Tuple{Int64,Int64}};
    node_vector::Vector{Vector{Float64}},
    num_quad_pts::Int64=2,
)
    constrained_dofs = get_constrained_dofs(constrained_dofs, length(connections), size(node_vector,1), 2)
    return Dataproblem(0.0, area, force, length(connections), size(node_vector,1), alpha, connections, constrained_dofs, node_vector, num_quad_pts, 2)
end




function fixedBarproblem1D(
    length::Float64,
    area::Float64,
    force::Function,
    num_ele::Int64,
    alpha::Float64;
    node_vector::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    num_quad_pts::Int64=2,
    right_fixed::Bool=true,
)
    num_node = num_ele + 1

    if isempty(node_vector)
        node_vector = [[x] for x in LinRange(0.0, length, num_node)]
    end
    constraints = [(1, 1)]
    if right_fixed
        push!(constraints, (num_node, 1))
    end
    constrained_dofs = get_constrained_dofs(constraints, num_ele, num_node, 1)
    connections = [(i, i + 1) for i in 1:num_ele]

    return Dataproblem(length, area, force, num_ele, num_node, alpha, connections, constrained_dofs, node_vector, num_quad_pts, 1)
end


function get_ndofs(problem::Dataproblem)
    ndof_u = ndof_lambda = problem.num_node * problem.dims
    ndof_e = ndof_s = ndof_mu = problem.num_ele
    return [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]
end
