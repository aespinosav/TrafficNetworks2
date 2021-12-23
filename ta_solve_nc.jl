#=
Functions to build and solve STAP for multiple OD pairs.
the 'nc' in the file name refers to 'network copies' since we use
this formulation explicitly (see Patriksson 1996).
=#

"""
Generate and solve multi pair stap by making copies of
the network (general expression of stap)
"""
function multi_pair_stap_nc(rn, ods, demands; regime=:ue)
    C = (regime == :ue ? 0.5 : 1.0)

    n = nv(rn.g)
    m = ne(rn.g)
    # Number of OD pairs
    n_ods = length(ods)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])  

    d_vects = SparseVector{Float64,Int64}[]
    for i in 1:length(ods)
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]
        push!(d_vects, d_vec)
    end

    A = incidence_matrix(rn.g)

    stap = Model(Gurobi.Optimizer)
    #set_silent(stap)

    # OD specific link flows
    @variable(stap, x[1:m,1:n_ods] >= 0)
    # Aggregate link flows (for expressing objective)
    @variable(stap, link_flow[1:m] >= 0)

    @constraint(stap, [i in 1:n_ods], A*x[:,i] .== d_vects[i])
    # Link flows have to add up
    @constraint(stap,
                inter_var_con[i in 1:m], # name of constraint
                link_flow[i] == sum(x[i,j] for j in 1:n_ods))

    # Objective function
    #For some reason this stopped working...
    #@objective(stap, Min, dot(a,link_flow) + (C*link_flow'*B*link_flow))
    @objective(stap, Min, sum(a.*link_flow) + C*link_flow'*B*link_flow )
    optimize!(stap)
    value.(x)
end
