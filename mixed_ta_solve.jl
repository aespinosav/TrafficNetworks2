#=
Functions to build and solve STAP for multiple OD pairs.
the 'nc' in the file name refers to 'network copies' since we use
this formulation explicitly (see Patriksson 1996).
=#

"""
Generate and solve multi pair stap by making copies of
the network (general expression of stap)
"""
function mixed_stap(rn, ods, demands, γ; ϵ=10^6)

    n = nv(rn.g)
    m = ne(rn.g)
    # Number of OD pairs
    n_ods = length(ods)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])  

    d_vects_hv = SparseVector{Float64,Int64}[]
    d_vects_av = SparseVector{Float64,Int64}[]
    for i in 1:length(ods)
        s, d = ods[i]
        d_vec = spzeros(n)
        d_vec[s] = -demands[i]
        d_vec[d] = demands[i]

        push!(d_vects_hv, (1-γ).*d_vec)
        push!(d_vects_av, γ.*d_vec)
    end

    A = incidence_matrix(rn.g)

    stap = Model(Gurobi.Optimizer)
    #set_silent(stap)

    # OD specific link flows HV
    @variable(stap, x[1:m,1:n_ods] >= 0)
    # OD specific link flows HV
    @variable(stap, y[1:m,1:n_ods] >= 0)

    # Aggregate link flows (for expressing objective)
    #@variable(stap, link_flow_[1:m] >= 0)

    # Link flows have to add up
    #@constraint(stap,
    #            inter_var_con[i in 1:m], # name of constraint
    #            link_flow[i] == sum(x[i,j] for j in 1:n_ods)
    #            )

    @constraint(stap,
                conservation_hvs[i in 1:n_ods],
                A*x[:,i] .== d_vects_hv[i])
    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*x[:,i] .== d_vects_av[i])

   hv_flow = sum(x, dims=2)
   av_flow = sum(y, dims=2)
   agg_flow = av_flow + hv_flow

   objective_hv = sum(hv_flow[i]*(a[i] + B[i,i]*av_flow[i]) +
                      0.5B[i,i]*hv_flow[i]^2 for i in 1:m)

   objective_av = sum(agg_flow[i]*a[i] + B[i,i]*agg_flow[i]^2
                      for i in 1:n)

   @objective(stap,
              Min,
              objective_hv
              )
   LBD = 0
   error = Inf
   while error >= epsilon


   end

    # Objective function
    #@objective(stap, Min, dot(a,link_flow) + (C*link_flow'*B*link_flow))


    optimize!(stap)
    value.(x)
end
