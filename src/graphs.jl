using Graphs, GraphPlot
import PlotlyJS

"""
    vivi_graph(inputs::Vector{Resource},techs::Vector{Tech},outputs::Vector{Resource},utils::Vector{Tech},store::Vector{Storage})
    Plots a graph representation
"""
function vivi_graph(inputs::Vector{Resource},techs::Vector{Tech},outputs::Vector{Resource},utils::Vector{Tech},store::Vector{Storage})
    g = DiGraph()
    vertex = Dict{String,Int64}()
    names = []

    # Techs and utils
    techs_all=vcat(techs,utils) # join techs and utils
    count = 1
    for tech in techs_all
        add_vertex!(g)
        vertex[tech.n] = count
        push!(names,tech.n)
        count += 1
    end
    n_tech = nv(g)

    for s in store
        add_vertex!(g)
        vertex[s.n] = count
        push!(names,s.n)
        count += 1
    end
    n_store = nv(g)

    resources = Set{String}()
    # There is an option to ommit this
    for tech in techs_all
        for input in tech.i
            push!(resources,input.t.n)
        end
        for output in tech.o
            push!(resources,output.t.n)
        end
    end

    # Resources in and out
    for input in inputs
        push!(resources,input.t.n)
    end
    for output in outputs
        push!(resources,output.t.n)
    end

    for resource in resources
        add_vertex!(g)
        vertex[resource] = count
        push!(names,resource)
        count += 1
    end

    # Create edge ##############################################################
    for tech in techs_all
        for input in tech.i
            i = vertex[input.t.n]
            j = vertex[tech.n]
            add_edge!(g,i,j)
        end
        for output in tech.o
            i = vertex[tech.n]
            j = vertex[output.t.n]
            add_edge!(g,i,j)
        end
    end

    for s in store
        if issubset([s.t.n],resources)
            i = vertex[s.n]
            j = vertex[s.t.n]
            add_edge!(g,i,j)
        end
    end

    pos_x, pos_y = GraphPlot.spring_layout(g)

    # Create plot points
    edge_x = []
    edge_y = []

    for edge in edges(g)
        push!(edge_x, pos_x[src(edge)])
        push!(edge_x, pos_x[dst(edge)])
        push!(edge_y, pos_y[src(edge)])
        push!(edge_y, pos_y[dst(edge)])
        push!(edge_x, NaN)
        push!(edge_y, NaN)
    end

    # Create edges
    edges_trace = PlotlyJS.scatter(
        mode="lines",
        x=edge_x,
        y=edge_y,
        line=PlotlyJS.attr(
            width=0.5,
            color="#888"
        ),
    )

    # Create nodes
    nodes_trace = PlotlyJS.scatter(
        x=pos_x[1:n_tech],
        y=pos_y[1:n_tech],
        mode="markers+text",
        text = [name for name in names[1:n_tech]],
        marker=PlotlyJS.attr(
            size=10,
            symbol="square",
            color="red",
        ),
        name = "Process",
    )

    # Create nodes
    nodes_trace2 = PlotlyJS.scatter(
        x=pos_x[n_tech+1:n_store],
        y=pos_y[n_tech+1:n_store],
        mode="markers+text",
        text = [name for name in names[n_tech+1:n_store]],
        marker=PlotlyJS.attr(
            size=10,
            symbol="diamond",
            color="yellow",
        ),
        name = "Storage",
    )

    nodes_trace3 = PlotlyJS.scatter(
        x=pos_x[n_store+1:end],
        y=pos_y[n_store+1:end],
        mode="markers+text",
        text = [name for name in names[n_store+1:end]],
        marker=PlotlyJS.attr(
            size=10,
            symbol="circle",
            color="#90EE90",
        ),
        name="Resource",
    )

    # Create Plot
    return PlotlyJS.plot(
        [edges_trace, nodes_trace, nodes_trace2, nodes_trace3],
        PlotlyJS.Layout(
            hovermode="closest",
            titlefont_size=16,
            showlegend=false,
            showarrow=false,
            xaxis=PlotlyJS.attr(showgrid=false, zeroline=false, showticklabels=false),
            yaxis=PlotlyJS.attr(showgrid=false, zeroline=false, showticklabels=false)
        )
    )
end

"""
    vivi_sankey(inputs::Vector{Resource},techs::Vector{Tech},outputs::Vector{Resource},utils::Vector{Tech},store::Vector{Storage};time::Int64=1,valueIndex::Int64=0,heatExergy::Bool=false,showHeat::Bool=true)
    Plots a sankey diagram
"""
function vivi_sankey(inputs::Vector{Resource},techs::Vector{Tech},outputs::Vector{Resource},utils::Vector{Tech},store::Vector{Storage};time::Int64=1,valueIndex::Int64=0,heatExergy::Bool=false,showHeat::Bool=true)
    # Could be just a set
    vertex = Dict{String,Int64}()
    names = []

    # Techs and utils
    techs_all=vcat(techs,utils) # join techs and utils
    count = 1
    for tech in techs_all
        vertex[tech.n] = count
        push!(names,tech.n)
        count += 1
    end
    n_s = count

    for s in store
        vertex[s.n] = count
        push!(names,s.n)
        count += 1
    end
    n_s2 = count

    # Resources
    resources = Set{String}()
    for tech in techs_all
        for input in tech.i
            push!(resources,input.t.n)
        end
        for output in tech.o
            push!(resources,output.t.n)
        end
    end
    for input in inputs
        push!(resources,input.t.n)
    end
    for output in outputs
        push!(resources,output.t.n)
    end
    for resource in resources
        vertex[resource] = count
        push!(names,resource)
        count += 1
    end

    # An error should be displayed if a tech and resource share the same name

    # Heat
    vertex["Heat network"] = count
    push!(names,"Heat network")

    # Connections
    loads_lims_per_tech = [load_limits(t) for t in techs_all]
    connections = []

    # Techs
    for (τ,tech) in enumerate(techs_all)
        size = tech.f
        load = round(tech.f_t[time]/size,digits=3)
        loads = loads_lims_per_tech[τ]
        if round(tech.f_t[time],digits=3) != 0
            s = 1
            while !(loads[s] <= load <=loads[s+1])
                s+=1
            end
            for input in tech.i
                i = vertex[input.t.n]
                j = vertex[tech.n]
                
                a,b = get_linear_parameters(tech,input.pw)
                index = length(input.r) == 1 ? 1 : time # Amount or size?
                qnt = input.r[index]*(tech.f_t[time]*a[s]+b[s]*size)
                if valueIndex != 0
                    value = length(input.t.c[:,valueIndex]) == 1 ? input.t.c[1,valueIndex] : input.t.c[time,valueIndex]
                    qnt*=value
                end
                append!(connections,[[i,j,qnt]])
            end
            
            for output in tech.o
                i = vertex[tech.n]
                j = vertex[output.t.n]

                a,b = get_linear_parameters(tech,output.pw)
                index = length(output.r) == 1 ? 1 : time # Amount or size?
                qnt = output.r[index]*(tech.f_t[time]*a[s]+b[s]*size)
                if valueIndex != 0
                    value = length(output.t.c[:,valueIndex]) == 1 ? output.t.c[1,valueIndex] : output.t.c[time,valueIndex]
                    qnt*=value
                end
                
                append!(connections,[[i,j,qnt]])
            end
        end
    end

    for s in store
        i = vertex[s.t.n]
        j = vertex[s.n]
        qnt = time >1 ? s.f_t[time]-s.f_t[time-1] : s.f_t[time] - s.a

        if qnt > 0
            append!(connections,[[i,j,qnt]])
        else
            append!(connections,[[j,i,-qnt]])
        end
    end

    # Heat
    n_heat = length(connections)
    heatNetwork = false
    if showHeat
        for (τ,tech) in enumerate(techs_all)
            
            size = tech.f
            load = round(tech.f_t[time]/size,digits=3)
            loads = loads_lims_per_tech[τ]

            if round(tech.f_t[time],digits=3) != 0
                s = 1
                while !(loads[s] <= load <=loads[s+1])
                    s+=1
                end
                
                qnt_in = 0
                qnt_out = 0
                for heat in tech.h

                    a,b = get_linear_parameters(tech,heat.pw)
                    heat_flow = heat.q*(tech.f_t[time]*a[s]+b[s]*size)

                    if !heatExergy 
                        qnt = heat_flow
                    else
                        qnt = heat_flow*(1-298.15/(heat.Tt-heat.Ts)*log(heat.Tt/heat.Ts))
                    end

                    if heat.Tt > heat.Ts
                        qnt_in += qnt
                    else
                        qnt_out += qnt
                    end
                    heatNetwork = true
                end
                append!(connections,[[count,vertex[tech.n],qnt_in]])
                append!(connections,[[vertex[tech.n],count,qnt_out]])
                    
            end
        end
    end

    # Create plot points
    edge_x = []
    edge_y = []
    amount = []
    for c in connections
        push!(edge_x, c[1]-1)
        push!(edge_y, c[2]-1)
        push!(amount, c[3])
    end

    # Color slices
    color = [i<n_s ? "blue" : (i<n_s2 ? "purple" : "#90EE90") for i=1:count]

    if heatNetwork
        color[count] = "red"
    end

    # Connection color
    color_2 = [i<=n_heat ? "rgba(186,186,186,0.5)" : "rgba(255,129,127,0.5)" for i=1:length(connections)] 

    return PlotlyJS.plot(PlotlyJS.sankey(
        node = PlotlyJS.attr(
          pad = 15,
          thickness = 20,
          line = PlotlyJS.attr(color = "black", width = 0.5),
          label = names,
          color = color
        ),
        link = PlotlyJS.attr(
          source = edge_x,
          target = edge_y,
          value = amount,
          color = color_2,
      ))
    )
end

# Graph shortcuts
vivi_graph(tech::Tech) = vivi_graph(tech.i,[tech],tech.o,Vector{Tech}(),Vector{Storage}())
vivi_graph(problem::Problem) = vivi_graph(problem.i,problem.p,problem.o,problem.ut,problem.st)

# Sankey shortcuts
vivi_sankey(tech::Tech;time::Int64=1,valueIndex::Int64=0,heatExergy::Bool=false,showHeat::Bool=true) = vivi_sankey(tech.i,[tech],tech.o,Vector{Tech}(),Vector{Storage}(),time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)
vivi_sankey(problem::Problem;time::Int64=1,valueIndex::Int64=0,heatExergy::Bool=false,showHeat::Bool=true) = vivi_sankey(problem.i,problem.p,problem.o,problem.ut,problem.st,time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)