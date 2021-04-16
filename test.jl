using LightGraphs, Test, Gurobi, TrafficNetworks2, SkeletonCities2, LinearAlgebra, SparseArrays

#const GRB_ENV = Gurobi.Env()

A = [0 1 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0]
g = SimpleDiGraph(A)
e_att = Dict(:a => [2, 4, 1, 4, 2], :b => [4, 2, 1, 2, 4])
n_att = Dict(:pos => [0.0 0.0; 1.0 1.0; 1.0 -1.0; 2.0 0.0])

rn = RoadNetwork(g, e_att, n_att)

od_single = od_matrix_from_pair(rn, 1, 4)

many_pairs = [(1,4), (1, 3), (2, 4)]
od_multi = multi_od_matrix(rn, many_pairs)
demands = [1.0, 0.5, 2.0]

multi_pair_stap(rn, od_multi, demands)
multi_pair_stap(rn, od_single, [0.5])

sd = StapData(rn, many_pairs, demands, :ue, true)
multi_pair_stap!(sd)

demand_ranges = [demands 2demands]
sd2 = StapData(rn, many_pairs, demand_ranges, :ue, true)
multi_pair_stap!(sd2)

mg = αβ_network_meta(10, 0.1, 1.5)
rn2 = skel2rn(mg)
dr = [3.0 6.0]
sd3 = StapData(rn2, [(1,45)], dr, :ue, true)
multi_pair_stap!(sd3)

drs = [3.0, 3.0]
sd4 = StapData(rn2, [(1,50), (10, 41)], drs, :ue, true)
multi_pair_stap!(sd4)

sd5 = StapData(rn2, [(1,50), (10, 41)], [drs 5drs], :ue, true)
multi_pair_stap!(sd5)
