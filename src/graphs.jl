using Graphs, GraphPlot
import PlotlyJS

function vivi_graph(inputs,techs,outputs,utils,store)
    g = DiGraph()
    vertex = Dict{String,Int64}()
    names = []

    # Techs and utils
    techs_all=vcat(techs,utils) # join techs and utils
    count = 1
    for tech in techs_all
        add_vertex!(g)
        vertex[tech.type] = count
        push!(names,tech.type)
        count += 1
    end
    n_tech = nv(g)

    for s in store
        add_vertex!(g)
        vertex["$(s.type) storage"] = count
        push!(names,"$(s.type) storage")
        count += 1
    end
    n_store = nv(g)

    resources = Set{String}()
    # There is an option to ommit this
    for tech in techs_all
        for input in tech.in
            push!(resources,input.type)
        end
        for output in tech.out
            push!(resources,output.type)
        end
    end

    # Resources in and out
    for input in inputs
        push!(resources,input.type)
    end
    for output in outputs
        push!(resources,output.type)
    end

    for resource in resources
        add_vertex!(g)
        vertex[resource] = count
        push!(names,resource)
        count += 1
    end

    # Create edge ##############################################################
    for tech in techs_all
        for input in tech.in
            i = vertex[input.type]
            j = vertex[tech.type]
            add_edge!(g,i,j)
        end
        for output in tech.out
            i = vertex[tech.type]
            j = vertex[output.type]
            add_edge!(g,i,j)
        end
    end

    for s in store
        if issubset([s.type],resources)
            i = vertex["$(s.type) storage"]
            j = vertex[s.type]
            add_edge!(g,i,j)
        end
    end

    # Names
    #graphplot(g,names=names,edgelabel=edgelabel_dict,nodeshape=:rect,curves=false)

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
            #title="Network Graph made with Julia",
            titlefont_size=16,
            showlegend=false,
            showarrow=false,
            xaxis=PlotlyJS.attr(showgrid=false, zeroline=false, showticklabels=false),
            yaxis=PlotlyJS.attr(showgrid=false, zeroline=false, showticklabels=false)
        )
    )
end
function vivi_sankey(inputs,techs,outputs,utils,store;time=1,valueIndex=0,heatExergy=false,showHeat=true)
    # Could be just a set
    vertex = Dict{String,Int64}()
    names = []

    # Techs and utils
    techs_all=vcat(techs,utils) # join techs and utils
    count = 1
    for tech in techs_all
        vertex[tech.type] = count
        push!(names,tech.type)
        count += 1
    end
    n_s = count

    # Resources
    resources = Set{String}()
    for tech in techs_all
        for input in tech.in
            push!(resources,input.type)
        end
        for output in tech.out
            push!(resources,output.type)
        end
    end
    for input in inputs
        push!(resources,input.type)
    end
    for output in outputs
        push!(resources,output.type)
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
    max_segments_per_tech,loads_lims_per_tech = segments_loads_per_tech(techs_all)

    connections = []
    for (τ,tech) in enumerate(techs_all)
        size = maximum(tech.size)
        load = tech.size[time]/size
        loads = loads_lims_per_tech[τ]

        if tech.size[time] != 0
            s = 1
            while !(loads[s] <= load <=loads[s+1])
                s+=1
            end
            
            for input in tech.in
                i = vertex[input.type]
                j = vertex[tech.type]
                
                a,b = get_linear_parameters(loads,input.loadEffect,max_segments_per_tech[τ])
                index = length(input.amount) == 1 ? 1 : time # Amount or size?
                qnt = input.amount[index]*(tech.size[time]*a[s]+b[s]*size)
                if valueIndex != 0
                    value = length(input.value[valueIndex]) == 1 ? input.value[valueIndex][1] : input.value[valueIndex][time]
                    qnt*=value
                end
                append!(connections,[[i,j,qnt]])
            end
            for output in tech.out
                i = vertex[tech.type]
                j = vertex[output.type]

                a,b = get_linear_parameters(loads,output.loadEffect,max_segments_per_tech[τ])
                index = length(output.amount) == 1 ? 1 : time # Amount or size?
                qnt = output.amount[index]*(tech.size[time]*a[s]+b[s]*size)
                if valueIndex != 0
                    value = length(output.value[valueIndex]) == 1 ? output.value[valueIndex][1] : output.value[valueIndex][time]
                    qnt*=value
                end
                
                append!(connections,[[i,j,qnt]])
            end
        end
    end

    n_heat = length(connections)
    heatNetwork = false
    if showHeat
        for (τ,tech) in enumerate(techs_all)
            
            size = maximum(tech.size)
            load = tech.size[time]/size
            loads = loads_lims_per_tech[τ]

            if tech.size[time] != 0
                s = 1
                while !(loads[s] <= load <=loads[s+1])
                    s+=1
                end
                
                for heat in tech.heat
                    i = heat.Tt > heat.Ts ? count : vertex[tech.type]
                    j = heat.Tt > heat.Ts ? vertex[tech.type] : count

                    a,b = get_linear_parameters(loads,heat.loadEffect,max_segments_per_tech[τ])
                    heat_flow = heat.h*(tech.size[time]*a[s]+b[s]*size)

                    if !heatExergy 
                        qnt = heat_flow # have a switch for exergy
                    else
                        qnt = heat_flow*(1-298.15/(heat.Tt-heat.Ts)*log(heat.Tt/heat.Ts))
                    end
                    append!(connections,[[i,j,qnt]])
                    heatNetwork = true
                end
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
    color = [i<n_s ? "blue" : "#90EE90" for i=1:count]
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
