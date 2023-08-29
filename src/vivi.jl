# Packages
using JuMP              # Optimization framework
using HiGHS             # MIT solver
using Gurobi            # Comercial solver
using DelimitedFiles    # To write files

# Importing functions
include("cascade.jl")   # Heat cascade
include("visual.jl")    # Visualization

#= ###################################
# Code structures
=# ###################################

mutable struct Resource
    type::String                    # Unique name
    amount::Vector{Real}            # Quantity that Techs will consume/produce
    unit::String                    # Unit of amount
    value::Vector{Vector{Real}}     # Specific value at each time step
end

function Resource(type,amount;unit="",value=[[0]])
    return Resource(type,amount,unit,value)
end

Resource(type::String,amount::Real,unit::String,value::Vector{Int64}) = Resource(type,[amount],unit,[[i] for i in value])   # Legacy code compatibility
Resource(type::String,amount::Real,unit::String,value::Vector{Float64}) = Resource(type,[amount],unit,[[i] for i in value]) # Legacy code compatibility
Resource(type::String,amount::Real;unit="",value=[[0]]) = Resource(type,[amount],unit=unit,value=value)

mutable struct Tech
    type::String                # Unique name
    in::Vector{Any}             # Consumed resources
    out::Vector{Any}            # Produced resources
    heat::Vector{HeatStruct}    # All heat transfers
    limits::Vector{Real}        # Max-min sizes
    cost::Vector{Real}          # Specific cost of technology - OBS: this may be on different basis
    loads::Vector{Real}         # Partial load limints (min,max)
    rate::Vector{Real}          # Ramping rate limits
    size::Vector{Real}          # Size factors - answer
end

function Tech(type;in=[],out=[],heat=[],limits=[0,Inf],cost=zeros(3),loads=[0,1],rate=[-1,1],size=[1])
    return Tech(type,in,out,heat,limits,cost,loads,rate,size)
end

Tech(type,in,out,heat;limits=[0,Inf],cost=zeros(3),loads=[0,1],rate=[-1,1]) = Tech(type,in=in,out=out,heat=heat,limits=limits,cost=cost,loads=loads,rate=rate) # Legacy code compatibility
Tech(type,in,out,heat,limits,cost::Real,loads;rate=[-1,1]) = Tech(type,in=in,out=out,heat=heat,limits=limits,cost=[cost,0],loads=loads,rate=rate) # Legacy code compatibility

mutable struct Storage
    type::String             # Resource name
    amount::Real             # Initial value
    max::Real                # Maximum storage capacity
    cost::Vector{Real}       # Specific cost of storage - OBS: this may be on different basis
    rate::Vector{Real}       # Max charging/discharge rate relative to storage capacity
    size::Vector{Real}       # Size factors - answer
end

function Storage(type,max;amount=0,cost=zeros(3),rate=[-1,1],size=[1])
    return Storage(type,amount,max,cost,rate,size)
end

Storage(type,amount,max,cost::Real,rate) = Storage(type,max,amount=amount,rate=rate,cost=[cost]) # Legacy code compatibility

mutable struct Problem
    inputs::Vector{Any}
    processes::Vector{Any}
    outputs::Vector{Any}
    utilities::Vector{Any}
    storage::Vector{Any}
end

function Problem(inputs,processes,outputs;utilities=[],storage=[])
    return Problem(inputs,processes,outputs,utilities,storage)
end

Problem(inputs,processes,outputs,utilities) = Problem(inputs,processes,outputs,utilities=utilities) # Legacy code compatibility

# Heat
Heat(h,Ts,Tt;dtExtra=0) = HeatStruct(h,Ts,Tt,dtExtra)
HeatStruct(h,Ts,Tt) = HeatStruct(h,Ts,Tt,0) # Legacy code compatibility

#= ###################################
# Optimization problem
=# ###################################

function vivi(problem::Problem;valueIndex=1,print=true,solver="HiGHS",capex=false,dt=15)
    # Get info #################################################################
    # Check if every tech was the same cost index lengths
 
    # Join techs and utilities
    techs=vcat(problem.processes,problem.utilities) # join techs and utils

    # Join every heat transfer
    t_q = techs[1].heat                 # stores the HeatStructs
    lims = [length(techs[1].heat)]      # saves the position of the HeatStructs for each Tech
    for i=2:length(techs)
        t_q = vcat(t_q,techs[i].heat)
        append!(lims,lims[i-1]+length(techs[i].heat))
    end

    # CHECK IF THE INPUTS AND OUTPUTS HAVE THE SAME TIME SIZE
 
    # Heat cascade constraint ##################################################
    # Temperature intervals
    Qk,locHot,locCold=heatCascade(t_q,[],[],forLP=true,dt=dt)

    # Define problem
    if solver == "HiGHS"
        LP=Model(HiGHS.Optimizer)
    else
        LP=Model(Gurobi.Optimizer)
    end
    set_silent(LP)

    # Define variables
    TkN=length(Qk[:,1])
    timeN=length(problem.inputs[1].amount)
    gammaN=length(techs)

    @variable(LP,R[1:TkN-1,1:timeN]>=0)
    @variable(LP,gammas[1:gammaN,1:timeN]>=0)

    for i=1:gammaN
        if techs[i].limits[1] != Inf
            for t = 1:timeN
                @constraint(LP,techs[i].limits[1] <= gammas[i,t] <= techs[i].limits[2])
            end
        end
    end

    # Define energy balance constraint
    @constraint(LP,con[k=1:TkN,t=1:timeN],0==0)
    for t=1:timeN
        for k=1:TkN-1
            set_normalized_coefficient(con[k,t],R[k,t],1)
        end
        for k=2:TkN
            set_normalized_coefficient(con[k,t],R[k-1,t],-1)
        end
    end
    
    # Heat cascade
    for t=1:timeN
        for k=1:TkN
            i = 1
            for b=1:gammaN
                if i<=lims[b] # Extra safety
                    set_normalized_coefficient(con[k,t],gammas[b,t],-sum(Qk[k,i:lims[b]]))
                    i = lims[b]+1
                else
                    set_normalized_coefficient(con[k,t],gammas[b,t],0)
                end
            end
        end
    end
    # OBS: normalized has to be the complete sum

    # Resource balance constraint ##############################################
    # List of resources
    resources = Set{String}()
    for input in problem.inputs
        push!(resources,input.type)
    end
    for output in problem.outputs
        push!(resources,output.type)
    end
    for tech in techs
        for input in tech.in
            push!(resources,input.type)
        end
        for output in tech.out
            push!(resources,output.type)
        end
    end
    # Improvement note: this implementation seems ugly and inefficient
    
    resN = length(resources)
    storeN = length(problem.storage)
    
    @variable(LP,resIN[1:resN,1:timeN]>=0)      # Improvement note: could only add variables that have been included as inputs
    @variable(LP,resOUT[1:resN,1:timeN]>=0)     # Improvement note: could only add variables that have been included as outputs
    
    # Resource index (name => index on constraint mass[i,t])
    valIN = zeros(resN,timeN) # coeffs in objective function
    valOUT = zeros(resN,timeN) # coeffs in objective function
    r_index = Dict{String,Int64}() # store resource name => index

    # Storage index (index on constraint mass[i,t] => index on )
    s_index = Dict{Int64,Int64}() # store resource name => index
    s_i = 1
    for s in problem.storage
        r_i = 1
        for r in resources
            if r == s.type
                s_index[r_i] = s_i
                break
            end
            r_i += 1
        end
        s_i += 1
    end

    # Intial storage constraint
    @variable(LP,0<=store_level[i=1:storeN,t=1:timeN+1]<=problem.storage[i].max)
    for i = 1:storeN
        @constraint(LP,store_level[i,1] == problem.storage[i].amount)
    end

    # Balance constraint
    @constraint(LP,mass[i=1:resN,t=1:timeN],resIN[i,t]-resOUT[i,t]==0)
    i = 1 # resource counter
    for r in resources
        # Techs
        for j=1:gammaN
            for res in techs[j].in
                if res.type == r
                    for t = 1:timeN
                        if length(res.amount) == 1
                            set_normalized_coefficient(mass[i,t],gammas[j,t],-res.amount[1]) # Assuming techs will not change with time
                        else
                            set_normalized_coefficient(mass[i,t],gammas[j,t],-res.amount[t]) # Assuming techs will change with time
                        end
                    end
                    break
                end
            end
            for res in techs[j].out
                if res.type == r
                    for t = 1:timeN
                        if length(res.amount) == 1
                            set_normalized_coefficient(mass[i,t],gammas[j,t],res.amount[1])
                        else
                            set_normalized_coefficient(mass[i,t],gammas[j,t],res.amount[t])
                        end
                    end
                    break
                end
            end
        end

        # Inputs
        isInput = false
        for input in problem.inputs
            if input.type == r
                isInput = true
                if length(input.value[valueIndex]) == timeN
                    # If is defined for each timestep
                    for t = 1:timeN
                        valIN[i,t] = input.value[valueIndex][t]
                        if input.amount[t] != Inf
                            @constraint(LP,resIN[i,t] == input.amount[t])
                        end
                    end
                elseif mod(timeN,length(input.value[valueIndex])) == 0
                    # If it's a multiple
                    n = length(input.value[valueIndex])
                    period = Int(timeN/n)
                    for t = timeN:-period:1
                        for tt = (t-period+1):t
                            valIN[i,tt] = input.value[valueIndex][n]
                        end
                        if input.amount[n] != Inf
                            @constraint(LP,sum(resIN[i,x] for x in (t-period+1):t) == input.amount[n])
                        end
                        n -= 1
                    end
                end
                # An error if the lengths are very different (actually should be earlier in the code)
                break
            end
        end
        if !isInput
            for t = 1:timeN
                @constraint(LP,resIN[i,t] == 0)
            end
        end

        # Outputs
        isOutput = false
        for output in problem.outputs
            if output.type == r
                isOutput = true
                if length(output.value[valueIndex]) == timeN
                    for t = 1:timeN
                        valOUT[i,t] = output.value[valueIndex][t]
                        if output.amount[t] != Inf
                            @constraint(LP,resOUT[i,t] == output.amount[t])
                        end
                    end
                elseif mod(timeN,length(output.value[valueIndex])) == 0
                    n = length(output.value[valueIndex])
                    period = Int(timeN/n)
                    for t = timeN:-period:1
                        # The value should be averaged
                        for tt = (t-period+1):t
                            valOUT[i,tt] = output.value[valueIndex][n]
                        end
                        if output.amount[n] != Inf
                            @constraint(LP,sum(resOUT[i,x] for x in (t-period+1):t) == output.amount[n])
                        end
                        n -= 1
                    end
                end
                break
            end
        end
        if !isOutput
            for t=1:timeN
                @constraint(LP,resOUT[i,t] == 0)
            end
        end

        # Storage
        for s in problem.storage
            if s.type == r
                for t = 1:timeN
                    set_normalized_coefficient(mass[i,t],store_level[s_index[i],t],1)
                    set_normalized_coefficient(mass[i,t],store_level[s_index[i],t+1],-1)
                end
            end
        end

        r_index[r] = i # store info
        i += 1 # resource counter
    end

    # Technology size constraint ##############################################
    @variable(LP,gamma_size[1:gammaN]>=0) # size of each technology
    @constraint(LP,highload[i=1:gammaN,t=1:timeN],gamma_size[i]*techs[i].loads[2]>=gammas[i,t])
    @constraint(LP,lowload[i=1:gammaN,t=1:timeN],gamma_size[i]*techs[i].loads[1]<=gammas[i,t])

    # Use or not of a technology
    @variable(LP, gamma_use[1:gammaN],Bin)
    @constraint(LP, useProcess[i=1:gammaN], gamma_size[i] <= gamma_use[i]*1E6) # not the best way to solve this actually because sizes can be higher

    # Ramping
    for i=1:length(techs)
        if techs[i].rate[1] != -1
            for t = 1:timeN-1
                @constraint(LP,gamma_size[i]*techs[i].rate[1] <= gammas[i,t]-gammas[i,t+1])
            end
        end
        if techs[i].rate[2] != 1
            for t = 1:timeN-1
                @constraint(LP,gamma_size[i]*techs[i].rate[2] >= gammas[i,t]-gammas[i,t+1])
            end
        end
    end

    # Storage size constraint ##############################################
    @variable(LP,gamma_size2[1:storeN]>=0) # size of each technology
    @constraint(LP,size2[i=1:storeN,t=1:timeN],gamma_size2[i]>=store_level[i,t]) # It should be higher or equal to the highest load

    for i=1:length(problem.storage)
        if problem.storage[i].rate[1] != -1
            for t = 1:timeN-1
                @constraint(LP,gamma_size2[i]*problem.storage[i].rate[1] <= store_level[i,t]-store_level[i,t+1])
            end
        end
        if problem.storage[i].rate[2] != 1
            for t = 1:timeN-1
                @constraint(LP,gamma_size2[i]*problem.storage[i].rate[2] >= store_level[i,t]-store_level[i,t+1])
            end
        end
    end

    # Set objective function ###################################################
    # Define objective function
    if capex
        @objective(LP,Max,sum(sum(resOUT[i,t]*valOUT[i,t]-resIN[i,t]*valIN[i,t] for i=1:resN) for t=1:timeN)-(sum(gamma_size[i]*techs[i].cost[1]+techs[i].cost[2]*gamma_use[i] for i=1:gammaN)+sum(gamma_size2[i]*problem.storage[i].cost[1] for i =1:storeN))*timeN)
#         @objective(LP,Max,sum(sum(resOUT[i,t]*valOUT[i,t]-resIN[i,t]*valIN[i,t] for i=1:resN) for t=1:timeN)-(sum(gamma_size[i]*techs[i].cost[1] for i=1:gammaN)+sum(gamma_size2[i]*problem.storage[i].cost[1] for i =1:storeN))*timeN)
    else
        @objective(LP,Max,sum(sum(resOUT[i,t]*valOUT[i,t]-resIN[i,t]*valIN[i,t] for i=1:resN) for t=1:timeN))
    end
    # Return results ###########################################################
    # Solve problem
    optimize!(LP)
    print && println(solution_summary(LP))
    
    gamma_opt = [value(gammas[s,t]) < 1E-3 ? 0 : value(gammas[s,t]) for s=1:gammaN for t=1:timeN]
    gamma_opt = reshape(gamma_opt,(timeN,gammaN))

    bin_use = [round(value(gamma_use[s]),sigdigits=5) for s=1:gammaN]

    store_opt = [round(value(store_level[s,t]),sigdigits=5) for s=1:storeN for t=1:timeN+1]
    store_opt = reshape(store_opt,(timeN+1,storeN))
    
    for i=1:storeN
        println(" Storage $(problem.storage[i].type) = $(maximum(store_opt[:,i]))")
    end
    
    # Print table
    larger = length("Tech")
    for tech in techs
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

        for i = 1:length(techs)
            n = " " * techs[i].type
            while length(n) < larger + 3
                n = n * " "
            end
            v = maximum(gamma_opt[:,i])
            println("$n | $v")
        end
    end

    # Return info ##############################################################
    inputs_ans = deepcopy(problem.inputs)
    outputs_ans = deepcopy(problem.outputs)
    processes_ans = deepcopy(problem.processes)
    utils_ans = deepcopy(problem.utilities)
    store_ans = deepcopy(problem.storage)

    # Inputs
    for input in inputs_ans
        input.amount=zeros(timeN)
        for t = 1:timeN
            ans_in = value(resIN[r_index[input.type],t])
            input.amount[t] = ans_in 
        end
    end

    # Output
    for output in outputs_ans
        output.amount=zeros(timeN)        
        for t = 1:timeN
            ans_out = value(resOUT[r_index[output.type],t])
            output.amount[t] = ans_out
        end
    end

    # Techs
    for i=1:length(processes_ans)
        processes_ans[i].size = gamma_opt[:,i]
    end

    # Utilities
    for i=(length(problem.processes)+1):length(gamma_opt[1,:])
        utils_ans[i-length(problem.processes)].size = gamma_opt[:,i]
    end

    # Storage
    for i=1:length(store_opt[1,:])
        store_ans[i].size = store_opt[:,i]
    end

    return Problem(inputs_ans,processes_ans,outputs_ans,utils_ans,store_ans)
end

# Graph shortcuts
vivi_graph(tech::Tech) = vivi_graph(tech.in,[tech],tech.out,[],[])
vivi_graph(problem::Problem) = vivi_graph(problem.inputs,problem.processes,problem.outputs,problem.utilities,problem.storage)

# Sankey shortcuts
vivi_sankey(tech::Tech;time=1,valueIndex=0,heatExergy=false,showHeat=true) = vivi_sankey(tech.in,[tech],tech.out,[],[],time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)
vivi_sankey(problem::Problem;time=1,valueIndex=0,heatExergy=false,showHeat=true) = vivi_sankey(problem.inputs,problem.processes,problem.outputs,problem.utilities,problem.storage,time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)

# Graphs
vivi_cc(problem::Problem) = vivi_plot(problem.processes,problem.utilities,"cc")
vivi_cc(techs::Tech) = vivi_plot([techs],[],"cc")
vivi_gcc(techs::Tech) = vivi_plot([techs],[],"gcc")

#=
vivi!(problem::Problem;valueIndex=1,Max=true,print=true)::Problem = vivi(problem;valueIndex=valueIndex,Max=Max,print=print,overwrite=true)
# Shortcut to plots
# Shortcuts
vivi_cc(t::Tech) = vivi_plot([t],[],"cc")
vivi_gcc(t::Tech) = vivi_plot([t],[],"gcc")
vivi_icc(p::Problem) = vivi_plot(p.processes,p.utilities,"icc")
vivi_gcc(p::Problem) = vivi_plot(p.processes,p.utilities,"gcc")
vivi_cc(p::Problem) = vivi_plot(p.processes,p.utilities,"cc")
=#