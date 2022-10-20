using Graphs, GraphRecipes

function vivi_graph(inputs,techs,outputs,utils)
    # Create vertex ############################################################
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

    # Resources in and out
    resources = Set{String}()
    for input in inputs
        push!(resources,input.type)
    end
    for output in outputs
        push!(resources,output.type)
    end
    #=
    for tech in techs_all
        for input in tech.in
            push!(resources,input.type)
        end
        for output in tech.out
            push!(resources,output.type)
        end
    end
    =#

    for resource in resources
        add_vertex!(g)
        vertex[resource] = count
        push!(names,resource)
        count += 1
    end

    # Create edge ##############################################################
    edgelabel_dict = Dict()
    for tech in techs_all
        for input in tech.in
            if issubset([input.type],resources)
                i = vertex[input.type]
                j = vertex[tech.type]
                add_edge!(g,i,j)
                edgelabel_dict[(i, j)] = "$(round(input.amount,sigdigits=3)) $(input.unit)"
            else
                for tech2 in techs_all
                    if tech.type != tech2.type
                        for output in tech2.out
                            if output.type == input.type
                                i = vertex[tech2.type]
                                j = vertex[tech.type]
                                add_edge!(g,i,j)
                                edgelabel_dict[(i, j)] = "$(round(input.amount,sigdigits=3)) $(input.unit)"
                                break
                            end
                        end
                    end
                end
            end

        end
        for output in tech.out
            if issubset([output.type],resources)
                i = vertex[tech.type]
                j = vertex[output.type]
                add_edge!(g,i,j)
                edgelabel_dict[(i, j)] = "$(round(output.amount,sigdigits=3)) $(output.unit)"

            else
                for tech2 in techs_all
                    if tech.type != tech2.type
                        for input in tech2.in
                            if input.type == output.type
                                i = vertex[tech.type]
                                j = vertex[tech2.type]
                                add_edge!(g,i,j)
                                edgelabel_dict[(i, j)] = "$(round(output.amount,sigdigits=3)) $(output.unit)"
                                break
                            end
                        end
                    end
                end
            end

        end
    end

    # Names
    return graphplot(g,names=names,edgelabel=edgelabel_dict,nodeshape=:rect,curves=false)
end
