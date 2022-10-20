# Description of the code

using JuMP
using GLPK

include("cascade.jl")
include("visual.jl")

"Mass resources e.g. fuel, water, power, etc."
mutable struct Resource
    type::String        # Unique name
    amount::Real        # Quantity that Techs will consume/produce
    unit::String          # Unit of amount (just for printing)
    value::Vector{Real} # Resource specific value [techno, economic, environment]
end
# OBS: Let the user handle the units (easier to implement)

"Conversion of resources"
struct Tech
    type::String                # Unique name
    in::Vector{Any}             # Consumed resources
    out::Vector{Any}            # Produced resources
    heat::Vector{HeatStruct}    # All heat transfers
end
# OBS : Power is handle as a resource (easier to implement)

"Problem structure"
struct Problem
    inputs::Vector{Any}
    processes::Vector{Any}
    outputs::Vector{Any}
    utilities::Vector{Any}
end

"Optimization problem formulation and solution"
function vivi(problem::Problem;valueIndex=1,Max=true,print=true,overwrite=false)::Problem
    # Get info #################################################################
    inputs = problem.inputs
    processes = problem.processes
    outputs = problem.outputs
    utils = problem.utilities

    # Join info ################################################################
    # Join techs and utilities
    techsN = length(processes)
    techs_all=vcat(processes,utils) # join techs and utils

    # Join every heat transfer
    t_q = techs_all[1].heat
    lims = [length(techs_all[1].heat)]
    for i=2:length(techs_all)
        t_q = vcat(t_q,techs_all[i].heat)
        append!(lims,lims[i-1]+length(techs_all[i].heat))
    end

    # Heat cascade constraint ##################################################
    # Temperature intervals
    Qk,locHot,locCold=heatCascade(t_q,[],[],forLP=true)

    # Define variables
    LP=Model(GLPK.Optimizer)
    TkN=length(Qk[:,1])
    @variable(LP,R[i=1:TkN]>=0)
    gammaN=length(techs_all)
    @variable(LP,gammas[i=1:gammaN]>=0)

    # Define energy balance constraint
    @constraint(LP,con[i=1:TkN],R[i]==0)
    for k=1:TkN
        k > 1 && set_normalized_coefficient(con[k],R[k-1],-1)
    end
    @constraint(LP,Rcon,R[TkN]==0)

    # Heat cascade
    for k=1:TkN
        i = 1
        for b=1:gammaN
            set_normalized_coefficient(con[k],gammas[b],-sum(Qk[k,i:lims[b]]))
            i = lims[b]+1
        end
    end
    # OBS: normalized has to be the complete sum

    # Resource balance constraint ##############################################
    # Variable
    resources = Set{String}()
    for input in inputs
        push!(resources,input.type)
    end
    for output in outputs
        push!(resources,output.type)
    end
    for tech in techs_all
        for input in tech.in
            push!(resources,input.type)
        end
        for output in tech.out
            push!(resources,output.type)
        end
    end
    resN = length(resources)
    @variable(LP,resIN[i=1:resN]>=0)
    @variable(LP,resOUT[i=1:resN]>=0)

    # Balance
    @constraint(LP,mass[i=1:resN],resIN[i]-resOUT[i]==0)
    i = 1 # resource counter
    valIN = zeros(resN) # coeffs in objective function
    valOUT = zeros(resN) # coeffs in objective function
    r_index = Dict{String,Int64}() # store resource name => index
    for r in resources
        # Techs
        for j=1:gammaN
            for res in techs_all[j].in
                if res.type == r
                    set_normalized_coefficient(mass[i],gammas[j],-res.amount)
                    break
                end
            end
            for res in techs_all[j].out
                if res.type == r
                    set_normalized_coefficient(mass[i],gammas[j],res.amount)
                    break
                end
            end
        end

        # Inputs
        isInput = false
        for input in inputs
            if input.type == r
                isInput = true
                valIN[i] = input.value[valueIndex]
                if input.amount != Inf
                    @constraint(LP,resIN[i] == input.amount)
                end
                break
            end
        end
        if !isInput
            @constraint(LP,resIN[i] == 0)
        end

        # Outputs
        isOutput = false
        for output in outputs
            if output.type == r
                isOutput = true
                valOUT[i] = output.value[valueIndex]
                if output.amount != Inf
                    @constraint(LP,resOUT[i] == output.amount)
                end
                break
            end
        end
        if !isOutput
            @constraint(LP,resOUT[i] == 0)
        end

        r_index[r] = i # store info
        i += 1 # resource counter
    end

    # Set objective function ###################################################
    # Define objective function
    if Max
        @objective(LP,Max,sum(resOUT[i]*valOUT[i]-resIN[i]*valIN[i] for i=1:resN))
    else
        @objective(LP,Max,sum(resOUT[i]*valOUT[i]-resIN[i]*valIN[i] for i=1:resN))
    end

    # Return results ###########################################################
    # Solve problem
    optimize!(LP)
    print && println(solution_summary(LP))
    gamma_opt = [round(value(gammas[s]),digits=5) for s=1:gammaN]

    # Print table
    larger = length("Tech")
    for tech in techs_all
        larger = length(tech.type) > larger ? length(tech.type) : larger
    end
    if print
        head = " Tech"
        line = "===="
        while length(head) < larger + 3
            head = head * " "
            line = line * "="
        end
        println("$head | Size factor")
        println("$line================")

        for i = 1:length(techs_all)
            n = " " * techs_all[i].type
            while length(n) < larger + 3
                n = n * " "
            end
            v = gamma_opt[i]
            println("$n | $v")
        end
    end

    # Return info ##############################################################
    if overwrite
        inputs_ans = inputs
        outputs_ans = outputs
        processes_ans = processes
        utils_ans = utils
    else
        inputs_ans = deepcopy(inputs)
        outputs_ans = deepcopy(outputs)
        processes_ans = deepcopy(processes)
        utils_ans = deepcopy(utils)
    end

    # Inputs
    for input in inputs_ans
        ans_in = value(resIN[r_index[input.type]])
        input.amount = ans_in >= 1E-7 ? ans_in : 0 # Avoid non-sense
    end

    # Output
    for output in outputs_ans
        ans_out = value(resOUT[r_index[output.type]])
        output.amount = ans_out >= 1E-7 ? ans_out : 0 # Avoid non-sense
    end

    # Techs
    for i=1:length(processes_ans)
        gamma = gamma_opt[i] > 1E-7 ? gamma_opt[i] : 0 # Avoid non-sense
        for input in processes_ans[i].in
            input.amount *= gamma
        end
        for output in processes_ans[i].out
            output.amount *= gamma
        end
        for heat in processes_ans[i].heat
            heat.h *= gamma
        end
    end

    # Utilities
    for i=(length(processes)+1):length(gamma_opt)
        gamma = gamma_opt[i] > 1E-7 ? gamma_opt[i] : 0 # Avoid non-sense
        for input in utils_ans[i-length(processes)].in
            input.amount *= gamma
        end
        for output in utils_ans[i-length(processes)].out
            output.amount *= gamma
        end
        for heat in utils_ans[i-length(processes)].heat
            heat.h *= gamma
        end
    end

    return Problem(inputs_ans,processes_ans,outputs_ans,utils_ans)
end

vivi!(problem::Problem;valueIndex=1,Max=true,print=true)::Problem = vivi(problem;valueIndex=valueIndex,Max=Max,print=print,overwrite=true)

# Shortcut to plots
vivi_graph(tech::Tech) = vivi_graph(tech.in,[tech],tech.out,[])
vivi_graph(problem::Problem) = vivi_graph(problem.inputs,problem.processes,problem.outputs,problem.utilities)
vivi_cc(t::Tech) = vivi_plot([t],[],"cc")
vivi_gcc(t::Tech) = vivi_plot([t],[],"gcc")
vivi_icc(p::Problem) = vivi_plot(p.processes,p.utilities,"icc")
vivi_gcc(p::Problem) = vivi_plot(p.processes,p.utilities,"gcc")
vivi_cc(p::Problem) = vivi_plot(p.processes,p.utilities,"cc")
