module TrafficNetworks2

using Graphs,
      LinearAlgebra,
      SparseArrays,
      JuMP,
      Gurobi,
      Ipopt,
      UnicodePlots,
      DataStructures

import Base.show,
       Graphs.incidence_matrix
       #LightGraphs.LinAlg.adjacency_matrix,
       #LightGraphs.LinAlg.incidence_matrix

export
    # From road_networks.jl
    RoadNetwork, random_od_pairs, od_matrix_from_pair, multi_od_matrix,
    adjacency_matrix, incidence_matrix,
    # From ta_solve_db.jl
    dest_nodes_flows, num_flows, make_od_mat_and_sort_d,
    make_demands_mat, multi_pair_stap,
    # From ta_solve_nc.jl
    multi_pair_stap_nc,
    # From stap_object.jl
    StapData, # multi_pair_stap!,
    # From mixed_ta_comb.jl
    mixed_stap_comb,
    # From edge_route_incidence.jl
    all_simple_paths, edge_route_incidence,
    # From plotting.jl
    term_plot

include("road_networks.jl")
include("ta_solve_db.jl")
include("ta_solve_nc.jl")
include("stap_object.jl")
include("mixed_ta_comb.jl")

# This file imports allsimplepaths.jl (not written by me)
include("edge_route_incidence.jl")
include("plotting.jl")

end
