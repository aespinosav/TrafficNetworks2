# Scratch pad for testing code

###
### Make OD pair and sort demands
###

###
### Multi pair (standard way)
###

###
### Mixed Equilibrium
###

m = Model(Gurobi.Optimizer)

@variable(m, x[1:3] >= 0)
@objective(m, Min, sum(x[i] for i in 1:3))

x_for_fixing = [2, 3, 4]
for i in 1:3
    fix(x[i], x_for_fixing[i], force=true)
end
optimize!(m)

for i in 1:length(x)
    unfix(x[i])
end
@constraint(m, x[1] >= 0)
optimize!(m)

fix(x[1], 2, force=true)
optimize!(m)
unfix(x[1])

equal_consts = all_constraints(stap, VariableRef, MOI.EqualTo{Float64})

#########################

grad_obj_x(u, v) = x -> sum((2*B[i,i]*(u[i]+v[i]) + a[i])*(x[i]-u[i]) for i in 1:length(x)) 


∇T_x(x, y) = [2*B[i,i]*(x[i]+y[i]) + a[i] for i in 1:length(x)]'
∇T_y(x, y) = [2*B[i,i]*(x[i]+y[i]) + 0.5*a[i] for i in 1:length(x)]'


#grad_obj_y = 

##############################################

adj_mat = [0 1 1 0;
           0 0 1 1;
           0 0 0 1;
           0 0 0 0]

a = [0.5, 1.0, 0.1, 1.0, 0.5]
b = [1.0, 0.5, 0.1, 0.5, 1.0]

g_directed = DiGraph(adj_mat)

positions = [0.05 0.5;
             0.5 0.75;
             0.5 0.25;
             0.95 0.5]

rn = RoadNetwork(g_directed,
                 Dict(:a=>a, :b=>b), 
                 Dict(:pos=>positions))

n = nv(rn.g)
m = ne(rn.g)

edge_route_incidence(rn, o, d)
