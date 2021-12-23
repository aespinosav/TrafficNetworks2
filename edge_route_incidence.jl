#using LightGraphs, SparseArrays
include("allsimplepaths.jl")

"""
Return array of simple paths (array of arrays), each path
is an array of integers, corresponding to the nodes
that are visited (in order).

    `all_paths_array(rn, o, d)`

    rn -- RoadNetwork
    o  -- origin (Int)
    d  -- destination (Int)
"""
function simple_paths_array(rn, o, d)
    simp_path_iter = all_simple_paths(rn.g, o, d)
    paths = collect(simp_path_iter)
end

"""
Returns the edge-route incidence matrix for a given 
road networ and set of paths
"""
function edge_route_incidence_matrix(rn, paths)

    sorted_edges = collect(edges(rn.g))
    m = ne(rn.g)
    r = length(paths)

    # Edge-Route incidence matrix
    Λ = spzeros(m, r)

    edge_index_paths = zeros.(Int, length.(paths) .- 1)
    for (i,p) in enumerate(paths)
        # Make Edge array to find indices
        e_p = [Edge(p[j],p[j+1]) for j in 1:length(p)-1]

        edge_indices = findfirst.(isequal.(e_p), [sorted_edges])
        edge_index_paths[i] .= edge_indices

        Λ[edge_indices, i] .= 1
    end

    Λ
end


"""
Get edge route incidence matrix and list of paths in the
indexed order they appear in the matrix.

Returns and E×R matrix, and an 2 arrays of paths.

The first path array are the paths as node sequences,
the second are as edge sequences.

Ordering of edges is what is given by collect(edges(g))
"""
function edge_route_incidence(rn, o, d)

    simp_path_iter = all_simple_paths(rn.g, o, d)
    paths = collect(simp_path_iter)

    sorted_edges = collect(edges(rn.g))
    m = ne(rn.g)
    r = length(paths)

    # Edge-Route incidence matrix
    Λ = spzeros(m, r)

    edge_index_paths = zeros.(Int, length.(paths) .- 1)
    for (i,p) in enumerate(paths)
        # Make Edge array to find indices
        e_p = [Edge(p[j],p[j+1]) for j in 1:length(p)-1]

        edge_indices = findfirst.(isequal.(e_p), [sorted_edges])
        edge_index_paths[i] .= edge_indices

        Λ[edge_indices, i] .= 1
    end

    Λ, edge_index_paths
end
