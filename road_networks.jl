"""
Road network object.

Contains a LightGraph simple directed graph (SimpleDiGraph),
and a dictionary with edge parameters.

Maybe not overparam the params
"""
mutable struct RoadNetwork
    g::SimpleDiGraph
    #edge_params::Dict{Symbol,Array{Float64,1}}
    #node_params::Dict{Symbol,Array{Float64,1}}
    edge_params
    node_params
end

#mutable struct RoadNetwork{T<:AbstractGraph}
#    g::T
#    #edge_params::Dict{Symbol,Array{Float64,1}}
#    #node_params::Dict{Symbol,Array{Float64,1}}
#    edge_params
#    node_params
#end

#function RoadNetwork(mg::MetaGraph)
#    g = SimpleDiGraph(mg)
#
#end

function show(io::IO, rn::RoadNetwork)
    output = "RoadNetwork:\nNodes - $(nv(rn.g)) \nEdges - $(ne(rn.g))\n"
    print(io, output)
end

adjacency_matrix(rn::RoadNetwork) = adjacency_matrix(rn.g)

incidence_matrix(rn::RoadNetwork) = incidence_matrix(rn.g)

"""
Returns an array of N OD pairs for a given graph g (uniformly chosen at random).
Nodes can be origins and destinations to multiple flows but not to the same flow.

    random_od_pairs(g::Graph, N)
"""
function random_od_pairs(g::AbstractGraph, N::Int)
    od_pair_array = zeros(Int, N, 2)
    for i in 1:N

        origin = rand(1:nv(g))
        destination = rand(1:nv(g))
        while destination == origin
            destination = rand(1:nv(g))
        end

        od_pair_array[i,:] = [origin, destination]
    end
    od_pair_array
end

random_od_pairs(rn::RoadNetwork, N::Int) = random_od_pairs(rn.g, N) 


"""
Generates OD matrix for graph 'g' for a single OD pair given as a 2-element tuple 'od_pair'.

Returns a sparse matrix.

    od_matrix_from_pair(g::AbstractGraph, od_pair::Tuple{Int64,Int64})
"""
function od_matrix_from_pair(g::AbstractGraph, o, d)
    n = nv(g)
    OD = spzeros(Int64, n, n)
    OD[o, d] = 1
    OD
end
od_matrix_from_pair(rn::RoadNetwork, o, d) = od_matrix_from_pair(rn.g, o, d)
od_matrix_from_pair(x::T, od::Array{Int,1}) where {T} = od_matrix_from_pair(x, od...)

"""
From a list of OD pairs, create OD matrix for multiple pairs
"""
function multi_od_matrix(g::AbstractGraph, ods)
    n = nv(g)
    OD = spzeros(Int64, n, n)

    for pair in ods
        o, d = pair
        OD += od_matrix_from_pair(g, o, d)
    end
    OD
end
multi_od_matrix(rn::RoadNetwork, ods) = multi_od_matrix(rn.g, ods)


