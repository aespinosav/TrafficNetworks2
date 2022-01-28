# Plotting networks to terminal

#=
using UnicodePlots
=#

function term_plot(rn::RoadNetwork)
    node_positions = rn.node_params[:pos]
    edge_num = ne(rn.g)

    edge_coords = Array{NTuple{4,Float64},1}(undef, edge_num)
    for (i, ed) in enumerate(edges(rn.g))
        x_src, y_src = node_positions[ed.src,:]
        x_dst, y_dst = node_positions[ed.dst,:]
        edge_coords[i] = (x_src, y_src, x_dst, y_dst)
    end

    canvas = BrailleCanvas(60,
                           25,
                           origin_x = 0.0,
                           origin_y = 0.0,
                           width = 1.0,
                           height = 1.0)
    for i in 1:edge_num
        lines!(canvas, edge_coords[i]...)
    end
    points!(canvas, node_positions[:,1], node_positions[:,2])
    canvas
end
