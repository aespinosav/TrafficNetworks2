"""
Struct for keeping tabs on stap solutions

demand_ranges should be a 2D array of nf×d_steps
"""
mutable struct StapData
    rn
    od_matrix
    od_pairs
    num_flow_vars
    d_steps
    demand_ranges
    regime
    dest_based
    flows
    agg_flows
end

function StapData(rn,
                  od_list::Array{T,1},
                  demand_ranges,
                  regime,
                  dest_based) where {T}
                  
    sorted_indices = sortperm(od_list)
    sorted_od_list = od_list[sorted_indices]
    OD = multi_od_matrix(rn, sorted_od_list)
    
    if length(size(demand_ranges))>1
        d_steps = size(demand_ranges)[2]
    else
        d_steps = 1
    end
    
    # Order demand ranges according to index of ordered ODs
    sorted_ranges = demand_ranges[sorted_indices,:]
    
    if dest_based == true
        nf = num_flows(OD)
    else
        nf = length(od_list)
    end
    
    flows = zeros(Float64, ne(rn.g), nf, d_steps)
    agg_flows = zeros(Float64, ne(rn.g), d_steps)
     
    StapData(rn,
             OD,
             sorted_od_list,
             nf,
             d_steps,
             sorted_ranges,
             regime,
             dest_based,
             flows,
             agg_flows)         
end

"""
Calculates traffic assignemnt for a given StapData object
that has been constructed with all relevant information
for the STAP to be solved.

Note that it should be flexible and be able to handle a 
single demand value or a demand range.
"""
function multi_pair_stap!(sd::StapData)
    
    for i in 1:sd.d_steps
        demands = sd.demand_ranges[:,i]
    
        x = multi_pair_stap(sd.rn, sd.od_matrix, demands, regime=sd.regime)
        
        sd.flows[:,:,i] = x
        sd.agg_flows[:,i] = sum(x, dims=2)
    end
end
