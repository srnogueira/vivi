# Packages
using JuMP              # Optimization framework
using HiGHS             # MIT solver
using Gurobi

# Importing functions
include("cascade.jl")   # Heat cascade

#= ###################################
# Code structures
=# ###################################

# TODO:
# - turn these default values into constraints without getting redefinition warnings (when the code is re-imported)
# - change value and loadEffect into Matrix structures

DEF_MAXSIZE = 1E3
DEF_COST = [[[0,0],[DEF_MAXSIZE,0]]]
DEF_VALUE = [[0]]
DEF_LOADEFF = []
DEF_UNIT = ""
DEF_LIMITS = [0,DEF_MAXSIZE]
DEF_RATE = [-1,1]

"""
Resource(type,amount,unit,value,loadEffect)

represents resources in the optimization problem
"""
mutable struct Resource
    type::String                    # Unique name
    amount::Vector{Real}            # Quantity that Techs will consume/produce
    unit::String                    # Unit of amount
    value::Vector{Vector{Real}}     # Specific value at each time step
    loadEffect::Vector{Vector{Real}}
end

function Resource(type,amount;unit=DEF_UNIT,value=DEF_VALUE,loadEffect=DEF_LOADEFF)
    return Resource(type,amount,unit,value,loadEffect)
end
Resource(type::String,amount::Real;unit="",value=DEF_VALUE) = Resource(type,[amount],unit=unit,value=value)
Resource(type::String,amount::Real,unit::String,value::Vector{Vector{Float64}}) = Resource(type,[amount],unit,value,[]) # Legacy code compatibility
Resource(type::String,amount::Real,unit::String,value::Vector{Vector{Int64}}) = Resource(type,[amount],unit,value,[]) # Legacy code compatibility
Resource(type::String,amount::Real,unit::String,value::Vector{Int64}) = Resource(type,[amount],unit,[[i] for i in value],[])   # Legacy code compatibility
Resource(type::String,amount::Real,unit::String,value::Vector{Float64}) = Resource(type,[amount],unit,[[i] for i in value],[]) # Legacy code compatibility
Resource(type::String,amount::Vector{Float64},unit::String,value::Vector{Vector{Int64}}) = Resource(type,amount,unit,value,[]) # Legacy code compatibility
Resource(type::String,amount::Vector{Float64},unit::String,value::Vector{Vector{Float64}}) = Resource(type,amount,unit,value,[]) # Legacy code compatibility


"""
Tech(type,in,out,heat,limits,cost,loads,rate,size)

"""
mutable struct Tech
    type::String                # Unique name
    in::Vector{Any}             # Consumed resources
    out::Vector{Any}            # Produced resources
    heat::Vector{HeatStruct}    # All heat transfers
    limits::Vector{Real}        # Max-min sizes
    cost          # Specific cost of technology - OBS: this may be on different basis
    loads::Vector{Real}         # Partial load limints (min,max)
    rate::Vector{Real}          # Ramping rate limits
    size::Vector{Real}          # Size factors - answer
end

function Tech(type;in=[],out=[],heat=[],limits=DEF_LIMITS,cost=DEF_COST,loads=[0,1,true],rate=DEF_RATE,size=[1])
    return Tech(type,in,out,heat,limits,cost,loads,rate,size)
end
Tech(type,in,out,heat;limits=DEF_LIMITS,cost=DEF_COST,loads=[0,1,true],rate=DEF_RATE) = Tech(type,in=in,out=out,heat=heat,limits=limits,cost=cost,loads=loads,rate=rate) # Legacy code compatibility
Tech(type,in,out,heat,limits,cost::Real,loads;rate=DEF_RATE) = Tech(type,in=in,out=out,heat=heat,limits=limits,cost=[cost,0],loads=loads,rate=rate) # Legacy code compatibility

mutable struct Storage # May not work as intended
    type::String             # Resource name
    amount::Real             # Initial value
    limits::Vector{Real}                # Maximum storage capacity - limits
    cost       # Specific cost of storage - OBS: this may be on different basis
    rate::Vector{Real}       # Max charging/discharge rate relative to storage capacity
    size::Vector{Real}       # Size factors - answer
end

function Storage(type,amount;limits=DEF_LIMITS,cost=DEF_COST,rate=DEF_RATE,size=[1])
    return Storage(type,amount,limits,cost,rate,size)
end
Storage(type,amount,max,cost::Real,rate) = Storage(type,amount,limits=[0,max],rate=rate,cost=[cost]) # Legacy code compatibility

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
Heat(h,Ts,Tt;dtExtra=0,loadEffect=DEF_LOADEFF) = HeatStruct(h,Ts,Tt,dtExtra,loadEffect)
HeatStruct(h,Ts,Tt;dtExtra=0,loadEffect=DEF_LOADEFF) = HeatStruct(h,Ts,Tt,dtExtra,loadEffect) # Legacy code compatibility

"""
add_tech_sizes!(model,techs)

adds the size of techs continuous variable and basic constraints
"""
function add_tech_sizes!(model,techs)
    n_techs = length(techs)
    @variable(model,f[τ=1:n_techs]) 
    @constraint(model,[τ=1:n_techs],f[τ]>=techs[τ].limits[1]) # size constraint
    @constraint(model,[τ=1:n_techs],f[τ]<=techs[τ].limits[2]) # size constraint
end

"""
create_capex_exp!(model,techs;valueIndex=1)

cretes a JuMP expression with the piecewise linearization of CAPEX
"""
function create_capex_exp!(model,techs;valueIndex=1)
    # Check number of segments for each Tech
    n_techs = length(techs)
    n_segments = Vector{Int}(undef,n_techs)

    for τ in eachindex(techs)
        points = techs[τ].cost[valueIndex]
        n_segments[τ] = length(points)-1 # Two points for every segment
        if n_segments[τ] <= 0
            throw("tech $(techs[τ].type) does not have cost points")
        end
    end

    # Create filter matrix
    matrix = Matrix{Int}(undef,n_techs,maximum(n_segments))
    for τ in eachindex(techs)
        for s=1:maximum(n_segments)
            matrix[τ,s] = s <= n_segments[τ] ? 1 : 0
        end
    end

    # Create binary variables
    @variable(model,y_s[τ=1:n_techs,s=1:maximum(n_segments); matrix[τ,s] > 0],Bin)
    @constraint(model,[τ=1:n_techs],sum(y_s[τ,:])<=1)

    # Create continuous variable
    add_tech_sizes!(model,techs)
    f = model[:f]

    @variable(model,f_s[τ=1:n_techs,s=1:maximum(n_segments); matrix[τ,s] > 0]>=0)
    @constraint(model,[τ=1:n_techs],sum(f_s[τ,:])==f[τ])
    for τ in eachindex(techs)
        points = techs[τ].cost[valueIndex]
        for s=1:n_segments[τ]
            @constraint(model,points[s][1]*y_s[τ,s]<=f_s[τ,s])
            @constraint(model,f_s[τ,s]<=points[s+1][1]*y_s[τ,s])
        end
    end

    # OBS: here is the point where additional costs could be added on

    # Create capex expression
    a = zeros(n_techs,maximum(n_segments))
    b = zeros(n_techs,maximum(n_segments))
    for τ in eachindex(techs)
        points = techs[τ].cost[valueIndex]
        for s=1:n_segments[τ]
            a[τ,s] = (points[s+1][2]-points[s][2])/(points[s+1][1]-points[s][1])
            b[τ,s] = points[s+1][2]-a[τ,s]*points[s+1][1]
        end
    end

    @expression(model,capex,sum(f_s[τ,s]*a[τ,s]+y_s[τ,s]*b[τ,s] for (τ,s) in eachindex(y_s)))
end

"""
create_capex_simple!(model,techs;valueIndex=1)
"""
function create_capex_simple!(model,techs;valueIndex=1)
    # Check number of segments for each Tech
    n_techs = length(techs)

    # Create continuous variable
    add_tech_sizes!(model,techs)
    f = model[:f]

    @expression(model,capex,sum(f[τ]*techs[τ].cost[valueIndex][end][2]/techs[τ].cost[valueIndex][end][1] for τ= 1:n_techs))
end

"""
get_linear_parameters(loads,points,max_segments)

returns the slope and intersection parameters of each linear segment
"""
function get_linear_parameters(loads,points,max_segments)
    a = zeros(max_segments)
    b = zeros(max_segments)
    if !isempty(points) && !isempty(points[1])
        # Calculate segment coefficients from provided points
        a_points = zeros(length(points)-1)
        b_points = zeros(length(points)-1)
        for p = 1:length(points)-1
            a_points[p] = (points[p+1][2]-points[p][2])/(points[p+1][1]-points[p][1])
            b_points[p] = points[p+1][2]-a_points[p]*points[p+1][1]
        end

        # Match coefficients with the tech segments
        p = 1
        for s = 1:max_segments
            if loads[s] == points[p][1]
                a[s] = a_points[p]
                b[s] = b_points[p]
                p+=1
            else
                if loads[s] < points[p][1]
                    throw("undefined load point for heat in tech") # Improve this error
                end
                a[s] = a[s-1] 
                b[s] = b[s-1]
            end
        end
    else
        a .= 1
        b .= 0
    end
    return a,b
end

"""
segments_loads_per_tech(techs)

returns info about the segments from different resources and heats in a tech
"""
function segments_loads_per_tech(techs)
    n_techs=length(techs)
    max_segments_per_tech = Vector{Int}(undef,n_techs)
    loads_lims_per_tech = Vector{Any}(undef,n_techs)
    for (τ,tech) in enumerate(techs)
        loads = Set()

        for i in tech.in
            for point in i.loadEffect
                length(point)>0 && push!(loads,point[1])
            end
        end
        for o in tech.out
            for point in o.loadEffect
                length(point)>0 && push!(loads,point[1])
            end
        end
        for h in tech.heat
            for point in h.loadEffect
                length(point)>0 && push!(loads,point[1])
            end
        end

        if tech.loads[1] == tech.loads[2]
            max_segments_per_tech[τ] = 1
            loads_lims_per_tech[τ] = [tech.loads[1],tech.loads[2]]
        else
            push!(loads,tech.loads[1]) # min load
            push!(loads,tech.loads[2]) # max load

            max_segments_per_tech[τ] = length(loads)-1
            loads_lims_per_tech[τ] = sort(collect(loads)) # should the tech overwrite the resource?

            if loads_lims_per_tech[τ][begin] != tech.loads[1] || loads_lims_per_tech[τ][end] != tech.loads[2]
                throw("Error")
            end
        end
    end
    return  max_segments_per_tech,loads_lims_per_tech
end

"""
create_tech_reformulation!(model,techs,resources)

creates matrixes with the tech resources and heat expressions using the reformulation strategy
"""
function create_tech_reformulation!(model,techs,resources,n_time;dt=15)
    # Find number of segments for each resource and heat of each tech
    n_techs = length(techs)
    max_segments_per_tech, loads_lims_per_tech = segments_loads_per_tech(techs)
    max_segments = maximum(max_segments_per_tech)
    # OBS: at least one segment will be created

    # Create filter matrix
    matrix = Matrix{Int}(undef,n_techs,max_segments)
    for τ in eachindex(techs)
        for s=1:max_segments
            matrix[τ,s] = s <= max_segments_per_tech[τ] ? 1 : 0
        end
    end

    # Create the variables - Reformulation strategy (f_ts are defined outside)
    f = model[:f]
    f_t = model[:f_t]
    @variable(model,f_st[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0]>=0)
    @variable(model,Ψ[τ=1:n_techs,t=1:n_time]>=0)
    @variable(model,ξ[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0]>=0)
    @variable(model,δ[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],Bin)

    # Constraints - Reformulation strategy
    @constraint(model,ref1[τ=1:n_techs,t=1:n_time],Ψ[τ,t]==f[τ])
    @constraint(model,ref2[τ=1:n_techs,t=1:n_time],sum(δ[τ,:,t])<=1)

    for (τ,tech) in enumerate(techs)
        if tech.loads[3] == false
            @constraint(model,[t=1:n_time],sum(δ[τ,:,t])==1) # Can't be shut down
        end
    end
    
    @constraint(model,ref3[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],f_st[τ,s,t] <= ξ[τ,s,t]*loads_lims_per_tech[τ][s+1])
    @constraint(model,ref4[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],f_st[τ,s,t] >= ξ[τ,s,t]*loads_lims_per_tech[τ][s])

    @constraint(model,ref5[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],ξ[τ,s,t] <= δ[τ,s,t]*techs[τ].limits[2]) # It seems to be the total limits and not the segments
    @constraint(model,ref6[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],ξ[τ,s,t] >= δ[τ,s,t]*techs[τ].limits[1])

    @constraint(model,ref7[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],Ψ[τ,t] <= (1-δ[τ,s,t])*techs[τ].limits[2]+ξ[τ,s,t])
    @constraint(model,ref8[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],Ψ[τ,t] >= (1-δ[τ,s,t])*techs[τ].limits[1]+ξ[τ,s,t])

    @constraint(model,ref9[τ=1:n_techs,t=1:n_time],f_t[τ,t]==sum(f_st[τ,:,t]))

    # Mass contributions
    n_resources = length(resources)
    
    r_ts = Matrix{Any}(undef,n_resources,n_time)
    for i=1:n_resources
        for t=1:n_time
            r_ts[i,t] = @expression(model,0)
        end
    end

    # Probably more efficient to loop for each Tech, also should not repeat itself
    r = 0 # A bit dangerous
    for resource in resources
        r += 1
        for τ in eachindex(techs)
            tech = techs[τ]

            for i in eachindex(tech.in)
                inlet = tech.in[i]
                if inlet.type == resource
                    # Finding piewise representation
                    loads = loads_lims_per_tech[τ]
                    points = inlet.loadEffect
                    a,b = get_linear_parameters(loads,points,max_segments_per_tech[τ])

                    # Add resource balance
                    for t = 1:n_time
                        if length(inlet.amount) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.amount[1]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.amount[t]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        end
                    end
                end
            end

            for o in eachindex(tech.out)
                outlet = tech.out[o]
                if outlet.type == resource
                    # Finding piewise representation
                    loads = loads_lims_per_tech[τ]
                    points = outlet.loadEffect
                    a,b = get_linear_parameters(loads,points,max_segments_per_tech[τ])
                    
                    # Add resource balance
                    for t = 1:n_time
                        if length(outlet.amount) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.amount[1]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.amount[t]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        end
                    end
                end
            end
        end
    end

    # Heats
    heats = []
    for tech in techs
        heats = vcat(heats,tech.heat)
    end
    if !isempty(heats)
        Qk,locHot,locCold=heatCascade(heats,[],[],forLP=true,dt=dt)
        n_Tk=length(Qk[:,1])
    else
        Qk = []
        n_Tk = 0
    end  

    # Create container
    q_kt = Matrix{Any}(undef,n_Tk,n_time)
    q_kt .= @expression(model,0)
 
    # Iterate
    start = 0
    for (τ,tech) in enumerate(techs)
        for (h,heat) in enumerate(tech.heat)
            Qki = Qk[:,h+start]
            loads = loads_lims_per_tech[τ]
            points = heat.loadEffect
            a,b = get_linear_parameters(loads,points,max_segments_per_tech[τ])

            for k=1:n_Tk
                for t=1:n_time
                    q_kt[k,t] = @expression(model,q_kt[k,t]-Qki[k]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                end
            end
        end
        start += length(tech.heat)
    end

    return r_ts,q_kt
end

"""
create_tech_simplified!(model,techs,resources)

creates matrixes with the tech resources and heat expressions using the simplified strategy
"""
function create_tech_simplified!(model,techs,resources,n_time;dt=15)
    # Find number of segments for each resource and heat of each tech
    n_techs = length(techs)

    # Create the variables - Reformulation strategy (f_ts are defined outside)
    f = model[:f]
    f_t = model[:f_t]
    
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] <= techs[τ].limits[2]) # It seems to be the total limits and not the segments
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] >= techs[τ].limits[1]) # It seems to be the total limits and not the segments
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] <= f[τ]*techs[τ].loads[2]) # It seems to be the total limits and not the segments
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] >= f[τ]*techs[τ].loads[1]) # It seems to be the total limits and not the segments


    # Mass contributions
    n_resources = length(resources)
    
    r_ts = Matrix{Any}(undef,n_resources,n_time)
    for i=1:n_resources
        for t=1:n_time
            r_ts[i,t] = @expression(model,0)
        end
    end

    # Probably more efficient to loop for each Tech, also should not repeat itself
    r = 0 # A bit dangerous
    for resource in resources
        r += 1
        for τ in eachindex(techs)
            tech = techs[τ]

            for i in eachindex(tech.in)
                inlet = tech.in[i]
                if inlet.type == resource

                    # Add resource balance
                    for t = 1:n_time
                        if length(inlet.amount) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.amount[1]*f_t[τ,t])
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.amount[t]*f_t[τ,t])
                        end
                    end
                end
            end

            for o in eachindex(tech.out)
                outlet = tech.out[o]
                if outlet.type == resource
                    
                    # Add resource balance
                    for t = 1:n_time
                        if length(outlet.amount) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.amount[1]*f_t[τ,t])
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.amount[t]*f_t[τ,t])
                        end
                    end
                end
            end
        end
    end

    # Heats
    heats = []
    for tech in techs
        heats = vcat(heats,tech.heat)
    end
    if !isempty(heats)
        Qk,locHot,locCold=heatCascade(heats,[],[],forLP=true,dt=dt)
        n_Tk=length(Qk[:,1])
    else
        Qk = []
        n_Tk = 0
    end  

    # Create container
    q_kt = Matrix{Any}(undef,n_Tk,n_time)
    q_kt .= @expression(model,0)
 
    # Iterate
    start = 0
    for (τ,tech) in enumerate(techs)
        for (h,heat) in enumerate(tech.heat)
            Qki = Qk[:,h+start]
            for k=1:n_Tk
                for t=1:n_time
                    q_kt[k,t] = @expression(model,q_kt[k,t]-Qki[k]*f_t[τ,t])
                end
            end
        end
        start += length(tech.heat)
    end

    return r_ts,q_kt
end

"""
add_ramping!(model,techs,n_time)

add ramping constraints for the techs in the optimization model
"""
function add_ramping!(model,techs,n_time)
    f = model[:f]
    f_t = model[:f_t]
    
    for i in eachindex(techs)
        @constraint(model,[t=1:n_time-1; techs[i].rate[1]!= DEF_RATE[1]],f[i]*techs[i].rate[1] <= f_t[i,t]-f_t[i,t+1])
        @constraint(model,[t=1:n_time-1; techs[i].rate[1]!= DEF_RATE[1]],f[i]*techs[i].rate[2] >= f_t[i,t]-f_t[i,t+1])
    end
end

"""
number_of_time_points(problem::Problem)

returns max. number of time points of problem
"""
function number_of_time_points(problem::Problem)
    n_time = 0
    for resource in vcat(problem.inputs,problem.outputs)
        n_time = max(n_time,length(resource.amount))
    end
    return n_time
end


"""
set_of_resources(problem::Problem)

returns the resources in a problem as a Set
"""
function list_of_resources(problem::Problem)
    resources = Set{String}()
    for r in vcat(problem.inputs,problem.outputs)
        push!(resources,r.type)
    end
    for tech in vcat(problem.processes,problem.utilities)
        for r in vcat(tech.in,tech.out)
            push!(resources,r.type)
        end
    end
    return collect(resources)
end

# Capacity factor
function CRF(;dis=0.08,years=20,hours=8760/2,f_COM=0.09,inf=0.02,f_CAPEX=1)
    i = (1+dis)/(1+inf)-1
    return i*(1+i)^years/((1+i)^years-1)/hours+f_COM*f_CAPEX/hours
end

#= ###################################
# Optimization problem
=# ###################################

"""
vivi(problem::Problem;valueIndex=1,print=true,solver=HiGHS.Optimizer,capex=false,dt=15)

Creates and solves the optimization problem
"""
function vivi(problem::Problem;valueIndex=1,print=true,solver="HiGHS",capex=false,dt=15,tmax=3*60,returnModel=false,CRF=CRF(),gap=false,LP=false)
    # Get info #################################################################
 
    # Join techs and utilities
    techs=vcat(problem.processes,problem.utilities)

    # Define problem
    if solver == "Gurobi"
        solver = Gurobi.Optimizer
    else
        solver = HiGHS.Optimizer
    end
    model = Model(solver)
    set_silent(model)
    set_time_limit_sec(model,tmax)
    if gap
        set_optimizer_attribute(model, "MIPGap", 0.01)
    end

    # Define variables
    n_time=number_of_time_points(problem)
    resources = list_of_resources(problem)

    n_res = length(resources)
    n_techs = length(techs)
    n_store = length(problem.storage)

    # Technology size/load variables ##############################################
    # TODO: 
    # - change the indexes from τ to the tech types

    if capex
        if LP
            create_capex_simple!(model,vcat(techs,problem.storage),valueIndex=valueIndex)
        else
            create_capex_exp!(model,vcat(techs,problem.storage),valueIndex=valueIndex)
        end
    else
        add_tech_sizes!(model,vcat(techs,problem.storage))
    end
    f = model[:f]
    @variable(model,f_t[τ=1:n_techs+n_store,t=1:n_time+1]>=0)

    # Reformulation strategy
    if LP
        r_ts,q_kt = create_tech_simplified!(model,techs,resources,n_time,dt=dt) # heat expressions are also created here
    else
        r_ts,q_kt = create_tech_reformulation!(model,techs,resources,n_time,dt=dt) # heat expressions are also created here
    end

    # Resources constraints #######################################################
    # TODO:
    # - change the indexes from i to the resources types
    # - avoid repetition

    @constraint(model,mass[i=1:n_res,t=1:n_time],r_ts[i,t]==0)

    binary_inOut = zeros(n_res,2)
    for (r,res) in enumerate(resources)
        for input in problem.inputs
            if input.type == res
                binary_inOut[r,1] = 1
                break
            end
        end
        for output in problem.outputs
            if output.type == res
                binary_inOut[r,2] = 1
                break
            end
        end
    end
    
    @variable(model,resIN[r=1:n_res,1:n_time;binary_inOut[r,1]==1]>=0) # Create sets indexed with the names of resources
    @variable(model,resOUT[r=1:n_res,1:n_time;binary_inOut[r,2]==1]>=0)
    
    for r=1:n_res
        if binary_inOut[r,1] == 1
            for t=1:n_time
                set_normalized_coefficient(mass[r,t],resIN[r,t],1)
            end
        end
        if binary_inOut[r,2] == 1
            for t=1:n_time
                set_normalized_coefficient(mass[r,t],resOUT[r,t],-1)
            end
        end
    end

    r_index = Dict{String,Int64}() # store resource name => index

    valIN = zeros(n_res,n_time)
    valOUT = zeros(n_res,n_time)

    for (i,r) in enumerate(resources)
        # Inputs
        for input in problem.inputs
            if input.type == r
                if length(input.value[valueIndex]) == n_time
                    # If is defined for each timestep
                    for t = 1:n_time
                        valIN[i,t] = input.value[valueIndex][t]
                        if input.amount[t] != Inf
                            @constraint(model,resIN[i,t] == input.amount[t])
                        end
                    end
                elseif mod(n_time,length(input.value[valueIndex])) == 0
                    # If it's a multiple
                    n = length(input.value[valueIndex])
                    period = Int(n_time/n)
                    for t = n_time:-period:1
                        for tt = (t-period+1):t
                            valIN[i,tt] = input.value[valueIndex][n]
                        end
                        if input.amount[n] != Inf
                            @constraint(model,sum(resIN[i,x] for x in (t-period+1):t) == input.amount[n])
                        end
                        n -= 1
                    end
                else
                    throw("mismatched lenghts of arrays")
                end
                break
            end
        end

        # Outputs
        for output in problem.outputs
            if output.type == r
                if length(output.value[valueIndex]) == n_time
                    for t = 1:n_time
                        valOUT[i,t] = output.value[valueIndex][t]
                        if output.amount[t] != Inf
                            @constraint(model,resOUT[i,t] == output.amount[t])
                        end
                    end
                elseif mod(n_time,length(output.value[valueIndex])) == 0
                    n = length(output.value[valueIndex])
                    period = Int(n_time/n)
                    for t = n_time:-period:1
                        # The value should be averaged
                        for tt = (t-period+1):t
                            valOUT[i,tt] = output.value[valueIndex][n]
                        end
                        if output.amount[n] != Inf
                            @constraint(model,sum(resOUT[i,x] for x in (t-period+1):t) == output.amount[n])
                        end
                        n -= 1
                    end
                else
                    throw("mismatched lenghts of arrays")
                end
                break
            end
        end

        r_index[r] = i
    end

    # Storage #########################################################
    # Intial storage constraint
    if n_store >= 1
        @constraint(model,[i=1:n_store,t=1:n_time+1],f_t[n_techs+i,t]<=problem.storage[i].limits[2])
        @constraint(model,[i=1:n_store,t=1:n_time+1],f_t[n_techs+i,t]>=problem.storage[i].limits[1])
        @constraint(model,[i=1:n_store],f_t[n_techs+i,1] == f[n_techs+i]*problem.storage[i].amount)
        @constraint(model,[i=1:n_store],f_t[n_techs+i,n_time+1] == f[n_techs+i]*problem.storage[i].amount)
        @constraint(model,size2[i=1:n_store,t=1:n_time+1],f[n_techs+i]>=f_t[n_techs+i,t]) # It should be higher or equal to the highest load
    end

    # Balance of storage
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

    for s in problem.storage
        for (i,r) in enumerate(resources)
            if s.type == r
                for t = 1:n_time
                    set_normalized_coefficient(mass[i,t],f_t[n_techs+s_index[i],t],1)
                    set_normalized_coefficient(mass[i,t],f_t[n_techs+s_index[i],t+1],-1)
                end
            end
        end
    end
        
    # Heat cascade constraint ##################################################
    n_Tk = isempty(q_kt) ? 0 : length(q_kt[:,1])

    # Define energy balance constraint
    @variable(model,R[1:n_Tk-1,1:n_time]>=0)
    @constraint(model,con[k=1:n_Tk,t=1:n_time],q_kt[k,t]==0)
    for t=1:n_time
        for k=1:n_Tk-1
            set_normalized_coefficient(con[k,t],R[k,t],1)
        end
        for k=2:n_Tk
            set_normalized_coefficient(con[k,t],R[k-1,t],-1)
        end
    end

    # Other constraints #######################################################
    # Ramping
    add_ramping!(model,vcat(techs,problem.storage),n_time)

    # Set objective function ###################################################
    # Define objective function
    cost_out = @expression(model,sum(sum(resOUT[i,t]*valOUT[i,t] for i=1:n_res if binary_inOut[i,2]==1) for t=1:n_time))
    cost_in = @expression(model,sum(sum(resIN[i,t]*valIN[i,t] for i=1:n_res if binary_inOut[i,1]==1) for t=1:n_time))
    if capex
        capex_exp = model[:capex]
        @objective(model,Max,cost_out-cost_in-capex_exp*n_time*CRF)
    else
        cost_out
        @objective(model,Max,cost_out-cost_in)
    end
    
    # Solve the problem
    optimize!(model)

    # Return results ###########################################################
    print && println(solution_summary(model))
    
    gamma_opt = [round(value(f_t[τ,t]),sigdigits=5) for τ=1:n_techs for t=1:n_time+1]
    gamma_opt = reshape(gamma_opt,(n_time+1,n_techs))

    store_opt = [round(value(f_t[s,t]),sigdigits=5) for s=n_techs+1:n_techs+n_store for t=1:n_time+1]
    store_opt = reshape(store_opt,(n_time+1,n_store))
    
    # Print table
    larger = length("Tech")
    for tech in vcat(techs,problem.storage)
        larger = length(tech.type) > larger ? length(tech.type) : larger
    end
    
    if print
        head = " Tech"
        line = "===="
        subline = "----"
        while length(head) < larger + 3
            head = head * " "
            line = line * "="
            subline = subline * "-"
        end
        println("$head | Size factor")
        println("$line================")

        for i in eachindex(techs)
            n = " " * techs[i].type
            while length(n) < larger + 3
                n = n * " "
            end
            v = maximum(gamma_opt[:,i])
            println("$n | $v")
        end

        println("$subline----------------")
        for s in eachindex(problem.storage)
            n = " " * problem.storage[s].type
            while length(n) < larger + 3
                n = n * " "
            end
            v = maximum(store_opt[:,s])
            println("$n | $v")
        end
    end

    # Returns the JuMP model instead
    if returnModel
        return model
    end

    # Return info ##############################################################
    inputs_ans = deepcopy(problem.inputs)
    outputs_ans = deepcopy(problem.outputs)
    processes_ans = deepcopy(problem.processes)
    utils_ans = deepcopy(problem.utilities)
    store_ans = deepcopy(problem.storage)

    # Inputs
    for input in inputs_ans
        input.amount=zeros(n_time)
        for t = 1:n_time
            ans_in = value(resIN[r_index[input.type],t])
            input.amount[t] = ans_in 
        end
    end

    # Output
    for output in outputs_ans
        output.amount=zeros(n_time)        
        for t = 1:n_time
            ans_out = value(resOUT[r_index[output.type],t])
            output.amount[t] = ans_out
        end
    end

    # Techs
    for (i,ans) in enumerate(processes_ans)
        ans.size = gamma_opt[:,i]
    end

    # Utilities
    for (i,ans) in enumerate(utils_ans)
        ans.size = gamma_opt[:,i+length(problem.processes)]
    end

    # Storage
    for (i,ans) in enumerate(store_ans)
        ans.size = store_opt[:,i]
    end

    return Problem(inputs_ans,processes_ans,outputs_ans,utils_ans,store_ans)
end

# Visualization
include("visual.jl")    # Visualization

# Graph shortcuts
vivi_graph(tech::Tech) = vivi_graph(tech.in,[tech],tech.out,[],[])
vivi_graph(problem::Problem) = vivi_graph(problem.inputs,problem.processes,problem.outputs,problem.utilities,problem.storage)

# Sankey shortcuts
vivi_sankey(tech::Tech;time=1,valueIndex=0,heatExergy=false,showHeat=true) = vivi_sankey(tech.in,[tech],tech.out,[],[],time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)
vivi_sankey(problem::Problem;time=1,valueIndex=0,heatExergy=false,showHeat=true) = vivi_sankey(problem.inputs,problem.processes,problem.outputs,problem.utilities,problem.storage,time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)

# Graphs
function vivi_plot(techs,utils,type;time=1)
    # Preparation
    t_q=[]
    max_segments_per_tech,loads_lims_per_tech = segments_loads_per_tech(techs)

    for (τ,tech) in enumerate(techs)
        aux = deepcopy(tech.heat)
        
        size = maximum(tech.size)
        load = tech.size[time]/size
        
        if size !=0 && load !=0
            loads = loads_lims_per_tech[τ]
            index = 1
            while !(loads[index] <= load <=loads[index+1])
                index+=1
            end

            for heat in aux
                a,b = get_linear_parameters(loads,heat.loadEffect,max_segments_per_tech[τ])
                heat.h*=tech.size[time]*a[index]+b[index]*size
            end
            append!(t_q,aux)
        end
    end
    techsTq = length(t_q)

    max_segments_per_tech,loads_lims_per_tech = segments_loads_per_tech(utils)
    for (τ,tech) in enumerate(utils)
        aux = deepcopy(tech.heat) # required

        size = maximum(tech.size)
        load = tech.size[time]/size
        if size !=0 && load !=0
            loads = loads_lims_per_tech[τ]
            index = 1
            while !(loads[index] <= load <=loads[index+1])
                index+=1
            end

            for heat in aux
                a,b = get_linear_parameters(loads,heat.loadEffect,max_segments_per_tech[τ])
                heat.h*=tech.size[time]*a[index]+b[index]*size
            end
            append!(t_q,aux)
        end
    end

    # Calcuations
    hotDummy = [HeatStruct(0.0,3000.0,2999.0)]
    coldDummy = [HeatStruct(0.0,1.0,2.0)]

    if type == "icc"
        merh,merc=mer(t_q[1:techsTq],hotDummy,coldDummy)
        plot = plot_icc(t_q,merh,merc,techsTq)
    elseif type == "cc"
        merh,merc=mer(t_q,hotDummy,coldDummy)
        plot=plot_cc(t_q,hotDummy,coldDummy,merh,merc)
    elseif type == "gcc"
        merh,merc=mer(t_q,hotDummy,coldDummy)
        plot=plot_gcc(t_q,hotDummy,coldDummy,merh,merc)
    end

    return plot
end

vivi_icc(techs,utils;time=1) = vivi_plot(techs,utils,"icc",time=time)
vivi_gcc(techs,utils;time=1) = vivi_plot(techs,utils,"gcc",time=time)
vivi_cc(techs,utils;time=1) = vivi_plot(techs,utils,"cc",time=time)

vivi_cc(problem::Problem;time=1) = vivi_plot(problem.processes,problem.utilities,"cc",time=time)
vivi_cc(techs::Tech;time=1) = vivi_plot([techs],[],"cc",time=time)
vivi_gcc(techs::Tech;time=1) = vivi_plot([techs],[],"gcc",time=time)