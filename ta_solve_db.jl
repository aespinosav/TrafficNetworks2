#= 
Functions to build and solve STAP for multiple OD pairs. 
Flows are considered based on destination only, in order to
reduce the number of variables per link.
=# 


"""
Get destinations of all OD pairs. od is given as a sparse array
"""
function dest_nodes_flows(od)
    I, J, V = findnz(od)
    num_pairs = length(V)
    sort(unique(J))
end


"""
Number of destination-based flows
"""
function num_flows(od)
    length(dest_nodes_flows(od))
end


"""
Returns dictionary indexed by destination, of the 
origin nodes of OD pairs with that given destination 

Tree is probably not the right term...
"""
function dest_tree(od::SparseMatrixCSC)

    I, J, V = findnz(od)
    destination_nodes = sort(unique(J))
    
    tree = Dict{Int64,Array{Int64,1}}()
    for (i,k) in enumerate(destination_nodes) 
        #indices in  dest array J that correspond to destination k
        indices = findall(x->x==k, J)
        #corresponding index in origin array I
        origins_to_k = I[indices] # These should be sorted already
        
        tree[k] = origins_to_k
    end
    tree
end

"""
Constructs dictionary of destination-based origins. Keys are destination index.

This function returns both the dictionary and the destination nodes,
to avoid having to iterate and sort the dictionary keys when wanting
to know the actual destination nodes of the graphs themselves.
"""
function dest_tree_and_flows(od::SparseMatrixCSC)

    I, J, V = findnz(od)
    destination_nodes = sort(unique(J))
    
    tree = Dict{Int64,Array{Int64,1}}()
    for (i,k) in enumerate(destination_nodes) 
        #indices in  dest array J that correspond to destination k
        indices = findall(x->x==k, J)
        #corresponding index in origin array I
        origins_to_k = I[indices]
        
        tree[k] = origins_to_k
    end
    tree, destination_nodes
end


"""
Retruns a demand matrix where each column has, for each flow 
(labeled by detination, in ascending index order), the source 
and sink values needed for the flow conservation constraints.

    `make_demands_mat(g, od, demands)`
    
    g:          Graph (also works with rn)
    od:         OD matrix (unweigted by actual demands)
    demands:    vector od demand flows
    
demands is a vector that contains demands of each od pair...
Care must be taken in terms of the order

Care must be taken with the sign of the demand, whether it is originating
or terminating.

Note that the d_lk term in the node loop to make the d_mat might have some
issues, the index in demands was originally demands[k] but it has now been
changed. I think it does the right thing now.
"""
function make_demands_mat(g::AbstractGraph, od::SparseMatrixCSC, demands::Array{T,1} where {T})    
    
    n = nv(g)
    
    flow_tree, dest_indices = dest_tree_and_flows(od)
    nf = length(dest_indices)
    
    d_mat = spzeros(n, nf)
    for (i,k) in enumerate(dest_indices)
        origins_k = flow_tree[k]
        suma_k = 0
        for (j,l) in enumerate(origins_k)
            # Note that in d_mat oringin indices match node indices
            d_lk = demands[i]
            d_mat[l,i] = -d_lk
            
            suma_k += d_lk
        end
        d_mat[k,i] = suma_k
    end
    
    d_mat
end
make_demands_mat(rn::RoadNetwork, od::SparseMatrixCSC, demands) = make_demands_mat(rn.g, od, demands)

"""
Makes demands matrix and returns d_mat as well as 
the sorted array of demands
"""
function make_od_mat_and_sort_d(g::AbstractGraph, od_list::Array{T,1} where {T}, demands)
    sorted_indices = sortperm(od_list)
    sorted_demands = demands[sorted_indices]
    od_mat = multi_od_matrix(g, od_list)
    
    od_mat, sorted_demands
end

"""
Solves STAP for multiple OD pairs. It does so by only keeping tabs on flows on a link according to the destination the flows have. So we lose information
on which origins things came from.
    
    ´multi_pair_stap(rn, od, demands; regime=:ue)´
"""
function multi_pair_stap(rn, od, demands; regime=:ue)
    C = (regime == :ue ? 0.5 : 1.0)
    
    n = nv(rn.g)
    m = ne(rn.g)
    
    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])
    
    d_mat = make_demands_mat(rn, od, demands)
    nf = size(d_mat)[2]
    
    stap_multi = Model(Gurobi.Optimizer)
    set_silent(stap_multi) # supress output
    
    # Disaggregate link flows (destination based)
    @variable(stap_multi, x[1:m,1:nf] >= 0)
    # Aggregate link variables
    @variable(stap_multi, link_flow[1:m] >= 0) 

    # Conservation constraint
    @constraint(stap_multi, incidence_matrix(rn.g)*x .== d_mat)
    # Link flows have to add up
    @constraint(stap_multi,
                inter_var_con[i in 1:m], # name of constraint
                link_flow[i] == sum(x[i,j] for j in 1:nf)
                )

    # Note how the link_flow variables are defined through the
    # Aggregate flow constratins

    # Objective function
    @objective(stap_multi, Min, dot(a,link_flow) + (C*link_flow'*B*link_flow))

    optimize!(stap_multi)
    
    #Maybe different things need to be returned
    value.(x)
end

#function multi_pair_stap(rn, od_list::Array{T,1} where {T}, demands; regime=:ue)
#    od_mat, sorted_demands = make_demands_mat_and_sort_d(rn.g, od_list, demands)
#    multi_pair_stap(rn, od_mat, sorted_demands, regime)
#end
