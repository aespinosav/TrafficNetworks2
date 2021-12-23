#=
In this file we implement in JuMP the mixed equilibrium model in which there are 2 user classes, one solving for UE and one for SO, as expressed by (van Vuren, van Vliet and Smith 1990) where it can all be summarised into a single cost function.
=#

# using LinearAlgebra, LightGraphs, SparseArrays, TrafficNetworks2, JuMP, Gurobi

function mixed_stap_comb_diag(rn,
                         ods,
                         demands,
                         γ;
                         tolerance=10^-6,
                         max_iters=50)

    # Number of OD pairs
    n_ods = length(ods)

    n = nv(rn.g)
    m = ne(rn.g)

    A = incidence_matrix(rn.g)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])

    me_obj_func(x, y) = sum( 0.5*B[i,i]*(x[i]+y[i])^2 +
                             a[i]*x[i] + y[i]*a[i]*0.5 for i in 1:length(x))

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

    #########################################

    stap = Model(Gurobi.Optimizer)
    set_silent(stap)

    # OD specific link flows HV
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, x .>= 0)
    # OD specific link flows HV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, y .>= 0)

    hv_flow = sum(x, dims=2)
    av_flow = sum(y, dims=2)
    #agg_flow = av_flow + hv_flow

    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])


    #mixed_objective = me_obj_func(x, y)
    mixed_objective = me_obj_func(hv_flow, av_flow)
    @objective(stap,
               Min,
               mixed_objective) 

    ### Initial values

    ue_sol = multi_pair_stap_nc(rn, ods, demands)

    # start with av's fixed
    y_fix_value = γ*ue_sol
    y_old = sum(y_fix_value, dims=2)[:]
    for i in 1:m
        for j in 1:n_ods
            fix(y[i,j], y_fix_value[i,j], force=true)
        end
    end

    x_start_val = (1-γ)*ue_sol
    x_old = sum(x_start_val, dims=2)[:]
    for i in 1:m
        for j in 1:n_ods
            set_start_value(x[i,j], x_start_val[i,j])
        end
    end

    ### Functions for term criterion

    ∇T_x(x, y) = [2*B[i,i]*(x[i]+y[i]) + a[i] for i in 1:length(x)]'
    ∇T_y(x, y) = [2*B[i,i]*(x[i]+y[i]) + 0.5*a[i] for i in 1:length(y)]'

    ### Optimisation

    obj_value_old = me_obj_func((1-γ)*ue_sol, γ*ue_sol)
    obj_diff = Inf
    n_iters = 0
    LBD = 0
    factors = Float64[]
    LB_check = Float64[]
    obj_array = Float64[]
    while obj_diff > tolerance && n_iters <= max_iters
    #while n_iters <= max_iters
        n_iters += 1

        optimize!(stap)

        ### Fix and unfix (Diagonalisation)

        for i in 1:m
            for j in 1:n_ods
                unfix(y[i,j]) 
            end
        end
        @constraint(stap, y .>= 0)

        x_val = value.(x)
        for i in 1:m
            for j in 1:n_ods
                fix(x[i,j], x_val[i,j], force=true)
            end
        end

        optimize!(stap)

        ### Fix and unfix

        for i in 1:m
            for j in 1:n_ods
                unfix(x[i,j])
            end
        end
        @constraint(stap, x .>= 0)

        y_val = value.(y)
        for i in 1:m
            for j in 1:n_ods
                fix(y[i,j], y_val[i,j], force=true)
            end
        end

        ### Check criteria

        obj_value_new = objective_value(stap)

        # Good criterion check
        x_new = sum(x_val, dims=2)[:]
        y_new = sum(y_val, dims=2)[:]

        T_lower = obj_value_old +
                  dot(∇T_x(x_old, y_old), (x_new - x_old)) +
                  dot(∇T_y(x_old, y_old), (y_new - y_old))

        #LBD = maximum([LBD, T_lower])

        conv_factor = (obj_value_new - T_lower)/T_lower
        push!(factors, conv_factor)
        push!(LB_check, obj_value_new > T_lower)
        push!(obj_array, obj_value_new)

        x_old = x_new
        y_old = y_new

        #Bad way of testing convergence
        obj_diff = obj_value_old - obj_value_new
        obj_value_old = obj_value_new
    end

    println("Total iterations: $(n_iters)")
    println("Last objective difference: $(obj_diff)")
    return value.(x), value.(y), factors, LB_check, obj_array
end



function mixed_stap_comb(rn,
                         ods,
                         demands,
                         γ;
                         tolerance=10^-6,
                         max_iters=50)

    # Number of OD pairs
    n_ods = length(ods)

    n = nv(rn.g)
    m = ne(rn.g)

    A = incidence_matrix(rn.g)

    a = rn.edge_params[:a]
    B = diagm(rn.edge_params[:b])

    me_obj_func(x, y) = sum( 0.5*B[i,i]*(x[i]+y[i])^2 +
                             a[i]*x[i] + y[i]*a[i]*0.5 for i in 1:length(x))

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

    #########################################

    stap = Model(Gurobi.Optimizer)
    set_silent(stap)

    # OD specific link flows HV
    @variable(stap, x[1:m,1:n_ods])
    @constraint(stap, x .>= 0)
    # OD specific link flows HV
    @variable(stap, y[1:m,1:n_ods])
    @constraint(stap, y .>= 0)

    hv_flow = sum(x, dims=2)
    av_flow = sum(y, dims=2)
    #agg_flow = av_flow + hv_flow

    @constraint(stap,
            conservation_hvs[i in 1:n_ods],
            A*x[:,i] .== d_vects_hv[i])

    @constraint(stap,
                conservation_avs[i in 1:n_ods],
                A*y[:,i] .== d_vects_av[i])


    #mixed_objective = me_obj_func(x, y)
    mixed_objective = me_obj_func(hv_flow, av_flow)
    @objective(stap,
               Min,
               mixed_objective) 

    optimize!(stap)

    return value.(x), value.(y)
end
