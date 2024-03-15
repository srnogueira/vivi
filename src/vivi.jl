# Packages
using JuMP,HiGHS,Gurobi

# Importing functions
include("cascade.jl")   # Heat cascade

#= ###################################
# Code structures
=# ###################################

MAXSIZE = 1E3

Base.@kwdef mutable struct ResourceType
    n::String = "-"
    c::Matrix{Real} = [1;;]
    u::String = "-"
end

"""
Resource(type,amount,unit,value,loadEffect)
    n: resource name
    r: rate
    u: unit
    c: cost  (t x type) [cost_1_t1 cost_2_t1; cost_1_t2 cost_2_t2]
    pw: piece-wise linearization

represents resources in the optimization problem
"""
Base.@kwdef mutable struct Resource
    t::ResourceType
    r::Vector{Real} = [0.0]
    pw::Matrix{Real} = [0 0; 1 1]
end

"""
Tech(type,in,out,heat,limits,cost,loads,rate,size)
    n: Name
    i: input
    o: output
    h: heat
    s: sizes
    c: costs (size - cost)
    l: loads
    rp: ramp
    f: size factor
    off : can be switch off?
"""
Base.@kwdef mutable struct Tech
    n::String = "-"
    i::Vector{Resource} = Array{Resource,1}()
    o::Vector{Resource} = Array{Resource,1}()
    h::Vector{Heat} = Array{Heat,1}()
    s::Vector{Real} = [0, MAXSIZE]
    c::Vector{Matrix{Real}} = [[0 0; MAXSIZE 0]]
    l::Vector{Real} = [0, 1]
    rp::Vector{Real} = [-1, 1]
    f::Vector{Real} = [1]
    off::Bool = true
end

"""
Tech(type,in,out,heat,limits,cost,loads,rate,size)
    n: Name
    a: amount
    h: heat
    s: sizes
    c: costs (size - cost)
    l: loads
    rp: ramp
    f: size factor
    
"""
Base.@kwdef mutable struct Storage
    n::String = "-"
    t::ResourceType = Resource()
    a::Real = 0
    s::Vector{Real} = [0, MAXSIZE]
    c::Vector{Matrix{Real}} = [[0 0; MAXSIZE 0]]
    rp::Vector{Real} = [-1, 1]
    f::Vector{Real} = [1]
end

Base.@kwdef mutable struct Problem
    i::Vector{Resource} = Array{Resource,1}()
    o::Vector{Resource} = Array{Resource,1}()
    p::Vector{Tech} = Array{Tech,1}()
    ut::Vector{Any} = Array{Tech,1}() 
    st::Vector{Any} = Array{Storage,1}()
end

"""
add_tech_sizes!(model,techs)

adds the size of techs continuous variable and basic constraints
"""
function add_tech_sizes!(model,techs)
    n_techs = length(techs)
    @variable(model,f[τ=1:n_techs]) 
    @constraint(model,[τ=1:n_techs],f[τ]>=techs[τ].s[1]) # size constraint
    @constraint(model,[τ=1:n_techs],f[τ]<=techs[τ].s[2]) # size constraint
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
        points = techs[τ].c[valueIndex]
        n_segments[τ] = length(points[:,1])-1 # Two points for every segment
        if n_segments[τ] <= 0
            throw("tech $(techs[τ].n) does not have cost points")
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
        points = techs[τ].c[valueIndex]
        for s=1:n_segments[τ]
            @constraint(model,points[s,1]*y_s[τ,s]<=f_s[τ,s])
            @constraint(model,f_s[τ,s]<=points[s+1,1]*y_s[τ,s])
        end
    end

    # OBS: here is the point where additional costs could be added on

    # Create capex expression
    a = zeros(n_techs,maximum(n_segments))
    b = zeros(n_techs,maximum(n_segments))
    for τ in eachindex(techs)
        points = techs[τ].c[valueIndex]
        for s=1:n_segments[τ]
            a[τ,s] = (points[s+1,2]-points[s,2])/(points[s+1,1]-points[s,1])
            b[τ,s] = points[s+1,2]-a[τ,s]*points[s+1,1]
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

    @expression(model,capex,sum(f[τ]*techs[τ].c[valueIndex][end,2]/techs[τ].c[valueIndex][end,1] for τ= 1:n_techs))
end

"""
get_linear_parameters(loads,points,max_segments)

returns the slope and intersection parameters of each linear segment
"""
function get_linear_parameters(loads,points,max_segments)
    a = zeros(max_segments)
    b = zeros(max_segments)
    if !isempty(points)
        # Calculate segment coefficients from provided points
        a_points = zeros(length(points[:,1])-1)
        b_points = zeros(length(points[:,1])-1)
        for p = 1:length(points[:,1])-1
            a_points[p] = (points[p+1,2]-points[p,2])/(points[p+1,1]-points[p,1])
            b_points[p] = points[p+1,2]-a_points[p]*points[p+1,1]
        end
        # Match coefficients with the tech segments
        p = 1
        for s = 1:max_segments
            if loads[s] == points[p,1] # Found a match
                a[s] = a_points[p]
                b[s] = b_points[p]
                p+=1
            else
                # (0) Is undefined
                if loads[s] < points[1,1]
                    # Should not happen
                    throw("Undefined load")
                elseif (points[p,1] <= loads[s] <= points[p+1,1])
                    # (A) is contained in a larger segment
                    a[s] =  a_points[p] 
                    b[s] = b_points[p]
                else
                    throw("Error")
                end

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
        for i in tech.i
            for point in i.pw[:,1]
                tech.l[1] <= point <= tech.l[2] && push!(loads,point)
            end
        end
        for o in tech.o
            for point in o.pw[:,1]
                tech.l[1] <= point <= tech.l[2] && push!(loads,point)
            end
        end
        for h in tech.h
            for point in h.pw[:,1]
                tech.l[1] <= point <= tech.l[2] && push!(loads,point)
            end
        end

        if tech.l[1] == tech.l[2]
            max_segments_per_tech[τ] = 1
            loads_lims_per_tech[τ] = [tech.l[1],tech.l[2]]
        else
            push!(loads,tech.l[1]) # min load
            push!(loads,tech.l[2]) # max load
            max_segments_per_tech[τ] = length(loads)-1
            loads_lims_per_tech[τ] = sort(collect(loads)) # should the tech overwrite the resource?
        end
    end
    return  max_segments_per_tech,loads_lims_per_tech
end

"""
create_tech_reformulation!(model,techs,resources)

creates matrixes with the tech resources and heat expressions using the reformulation strategy
"""
function create_tech_reformulation!(model,techs,resources,n_time)
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
        if tech.off == false
            @constraint(model,[t=1:n_time],sum(δ[τ,:,t])==1) # Can't be shut down
        end
    end
    
    @constraint(model,ref3[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],f_st[τ,s,t] <= ξ[τ,s,t]*loads_lims_per_tech[τ][s+1])
    @constraint(model,ref4[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],f_st[τ,s,t] >= ξ[τ,s,t]*loads_lims_per_tech[τ][s])

    @constraint(model,ref5[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],ξ[τ,s,t] <= δ[τ,s,t]*techs[τ].s[2]) # It seems to be the total limits and not the segments
    @constraint(model,ref6[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],ξ[τ,s,t] >= δ[τ,s,t]*techs[τ].s[1])

    @constraint(model,ref7[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],Ψ[τ,t] <= (1-δ[τ,s,t])*techs[τ].s[2]+ξ[τ,s,t])
    @constraint(model,ref8[τ=1:n_techs,s=1:max_segments,t=1:n_time;matrix[τ,s]>0],Ψ[τ,t] >= (1-δ[τ,s,t])*techs[τ].s[1]+ξ[τ,s,t])

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

            for i in eachindex(tech.i)
                inlet = tech.i[i]
                if inlet.t.n == resource
                    # Finding piewise representation
                    loads = loads_lims_per_tech[τ]
                    points = inlet.pw
                    a,b = get_linear_parameters(loads,points,max_segments_per_tech[τ])

                    # Add resource balance
                    for t = 1:n_time
                        if length(inlet.r) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.r[1]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.r[t]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        end
                    end
                end
            end

            for o in eachindex(tech.o)
                outlet = tech.o[o]
                if outlet.t.n == resource
                    # Finding piewise representation
                    loads = loads_lims_per_tech[τ]
                    points = outlet.pw
                    a,b = get_linear_parameters(loads,points,max_segments_per_tech[τ])
                    # Add resource balance
                    for t = 1:n_time
                        if length(outlet.r) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.r[1]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.r[t]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                        end
                    end
                end
            end
        end
    end

    # Heats
    heats = []
    for tech in techs
        heats = vcat(heats,tech.h)
    end
    if !isempty(heats)
        Qk=heatCascade(heats)
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
        for (h,heat) in enumerate(tech.h)
            Qki = Qk[:,h+start]
            loads = loads_lims_per_tech[τ]
            points = heat.pw
            a,b = get_linear_parameters(loads,points,max_segments_per_tech[τ])

            for k=1:n_Tk
                for t=1:n_time
                    q_kt[k,t] = @expression(model,q_kt[k,t]-Qki[k]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ]))
                end
            end
        end
        start += length(tech.h)
    end

    return r_ts,q_kt
end

"""
create_tech_simplified!(model,techs,resources)

creates matrixes with the tech resources and heat expressions using the simplified strategy
"""
function create_tech_simplified!(model,techs,resources,n_time)
    # Find number of segments for each resource and heat of each tech
    n_techs = length(techs)

    # Create the variables - Reformulation strategy (f_ts are defined outside)
    f = model[:f]
    f_t = model[:f_t]
    
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] <= techs[τ].s[2]) # It seems to be the total limits and not the segments
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] >= techs[τ].s[1]) # It seems to be the total limits and not the segments
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] <= f[τ]*techs[τ].l[2]) # It seems to be the total limits and not the segments
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] >= f[τ]*techs[τ].l[1]) # It seems to be the total limits and not the segments

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
            for i in eachindex(tech.i)
                inlet = tech.i[i]
                if inlet.t.n == resource
                    # Add resource balance
                    for t = 1:n_time
                        if length(inlet.r) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.r[1]*f_t[τ,t])
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]-inlet.r[t]*f_t[τ,t])
                        end
                    end
                end
            end

            for o in eachindex(tech.o)
                outlet = tech.o[o]
                if outlet.t.n == resource
                    
                    # Add resource balance
                    for t = 1:n_time
                        if length(outlet.r) == 1 # Maybe it does not make so much sense
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.r[1]*f_t[τ,t])
                        else
                            r_ts[r,t] = @expression(model,r_ts[r,t]+outlet.r[t]*f_t[τ,t])
                        end
                    end
                end
            end
        end
    end

    # Heats
    heats = []
    for tech in techs
        heats = vcat(heats,tech.h)
    end
    if !isempty(heats)
        Qk=heatCascade(heats)
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
        for (h,heat) in enumerate(tech.h)
            Qki = Qk[:,h+start]
            for k=1:n_Tk
                for t=1:n_time
                    q_kt[k,t] = @expression(model,q_kt[k,t]-Qki[k]*f_t[τ,t])
                end
            end
        end
        start += length(tech.h)
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
        @constraint(model,[t=1:n_time-1; techs[i].rp[1]!= -1],f[i]*techs[i].rp[1] <= f_t[i,t]-f_t[i,t+1])
        @constraint(model,[t=1:n_time-1; techs[i].rp[2]!= 1],f[i]*techs[i].rp[2] >= f_t[i,t]-f_t[i,t+1])
    end
end

"""
number_of_time_points(problem::Problem)

returns max. number of time points of problem
"""
function number_of_time_points(problem::Problem)
    n_time = 0
    for resource in vcat(problem.i,problem.o)
        n_time = max(n_time,length(resource.r))
    end
    return n_time
end


"""
set_of_resources(problem::Problem)

returns the resources in a problem as a Set
"""
function list_of_resources(problem::Problem)
    resources = Set{String}()
    for r in vcat(problem.i,problem.o)
        push!(resources,r.t.n)
    end
    for tech in vcat(problem.p,problem.ut)
        for r in vcat(tech.i,tech.o)
            push!(resources,r.t.n)
        end
    end
    return collect(resources)
end

#= ###################################
# Optimization problem
=# ###################################

"""
vivi(problem::Problem;valueIndex=1,print=true,solver=HiGHS.Optimizer,capex=false,dt=15)

Creates and solves the optimization problem
"""
function vivi(problem::Problem;valueIndex=1,print=true,solver="HiGHS",capex=false,tmax=3*60,returnModel=false,gap=false,LP=false)
    # Get info #################################################################
 
    # Join techs and utilities
    techs=vcat(problem.p,problem.ut)

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
    n_store = length(problem.st)

    # Technology size/load variables ##############################################
    # TODO: 
    # - change the indexes from τ to the tech types

    if capex
        if LP
            create_capex_simple!(model,vcat(techs,problem.st),valueIndex=valueIndex)
        else
            create_capex_exp!(model,vcat(techs,problem.st),valueIndex=valueIndex)
        end
    else
        add_tech_sizes!(model,vcat(techs,problem.st))
    end
    f = model[:f]
    @variable(model,f_t[τ=1:n_techs+n_store,t=1:n_time+1]>=0)

    # Reformulation strategy
    if LP
        r_ts,q_kt = create_tech_simplified!(model,techs,resources,n_time) # heat expressions are also created here
    else
        r_ts,q_kt = create_tech_reformulation!(model,techs,resources,n_time) # heat expressions are also created here
    end

    # Resources constraints #######################################################
    # TODO:
    # - change the indexes from i to the resources types
    # - avoid repetition

    @constraint(model,mass[i=1:n_res,t=1:n_time],r_ts[i,t]==0)

    binary_inOut = zeros(n_res,2)
    for (r,res) in enumerate(resources)
        for input in problem.i
            if input.t.n == res
                binary_inOut[r,1] = 1
                break
            end
        end
        for output in problem.o
            if output.t.n == res
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
        for input in problem.i
            if input.t.n == r
                if length(input.t.c[:,valueIndex]) == n_time
                    # If is defined for each timestep
                    for t = 1:n_time
                        valIN[i,t] = input.t.c[t,valueIndex]
                        if input.r[t] != Inf
                            @constraint(model,resIN[i,t] == input.r[t])
                        end
                    end
                else
                    throw("mismatched lenghts of arrays")
                end
                break
            end
        end

        # Outputs
        for output in problem.o
            if output.t.n == r
                if length(output.t.c[:,valueIndex]) == n_time
                    for t = 1:n_time
                        valOUT[i,t] = output.t.c[t,valueIndex]
                        if output.r[t] != Inf
                            @constraint(model,resOUT[i,t] == output.r[t])
                        end
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
        @constraint(model,[i=1:n_store,t=1:n_time+1],f_t[n_techs+i,t]<=problem.st[i].s[2])
        @constraint(model,[i=1:n_store,t=1:n_time+1],f_t[n_techs+i,t]>=problem.st[i].s[1])
        @constraint(model,[i=1:n_store],f_t[n_techs+i,1] == f[n_techs+i]*problem.st[i].a)
        @constraint(model,[i=1:n_store],f_t[n_techs+i,n_time+1] == f[n_techs+i]*problem.st[i].a)
        @constraint(model,size2[i=1:n_store,t=1:n_time+1],f[n_techs+i]>=f_t[n_techs+i,t]) # It should be higher or equal to the highest load
    end

    # Balance of storage
    s_index = Dict{Int64,Int64}() # store resource name => index
    
    s_i = 1
    for s in problem.st
        r_i = 1
        for r in resources
            if r == s.t.n
                s_index[r_i] = s_i
                break
            end
            r_i += 1
        end
        s_i += 1
    end

    for s in problem.st
        for (i,r) in enumerate(resources)
            if s.t.n == r
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
    add_ramping!(model,vcat(techs,problem.st),n_time)

    # Set objective function ###################################################
    # Define objective function
    cost_out = @expression(model,sum(sum(resOUT[i,t]*valOUT[i,t] for i=1:n_res if binary_inOut[i,2]==1) for t=1:n_time))
    cost_in = @expression(model,sum(sum(resIN[i,t]*valIN[i,t] for i=1:n_res if binary_inOut[i,1]==1) for t=1:n_time))
    if capex
        capex_exp = model[:capex]
        @objective(model,Max,cost_out-cost_in-capex_exp*n_time)
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
    for tech in vcat(techs,problem.st)
        larger = length(tech.n) > larger ? length(tech.n) : larger
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
            n = " " * techs[i].n
            while length(n) < larger + 3
                n = n * " "
            end
            v = maximum(gamma_opt[:,i])
            println("$n | $v")
        end

        println("$subline----------------")
        for s in eachindex(problem.st)
            n = " " * problem.st[s].t.n
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
    inputs_ans = deepcopy(problem.i)
    outputs_ans = deepcopy(problem.o)
    processes_ans = deepcopy(problem.p)
    utils_ans = deepcopy(problem.ut)
    store_ans = deepcopy(problem.st)

    # Inputs
    for input in inputs_ans
        input.r=zeros(n_time)
        for t = 1:n_time
            ans_in = value(resIN[r_index[input.t.n],t])
            input.r[t] = ans_in 
        end
    end

    # Output
    for output in outputs_ans
        output.r=zeros(n_time)        
        for t = 1:n_time
            ans_out = value(resOUT[r_index[output.t.n],t])
            output.r[t] = ans_out
        end
    end

    # Techs
    for (i,ans) in enumerate(processes_ans)
        ans.f = gamma_opt[:,i]
    end

    # Utilities
    for (i,ans) in enumerate(utils_ans)
        ans.f = gamma_opt[:,i+length(problem.p)]
    end

    # Storage
    for (i,ans) in enumerate(store_ans)
        ans.f = store_opt[:,i]
    end

    return Problem(i=inputs_ans,p=processes_ans,o=outputs_ans,ut=utils_ans,st=store_ans)
end

# Visualization
include("graphs.jl")

# Graph shortcuts
vivi_graph(tech::Tech) = vivi_graph(tech.i,[tech],tech.o,[],[])
vivi_graph(problem::Problem) = vivi_graph(problem.i,problem.p,problem.o,problem.ut,problem.st)

# Sankey shortcuts
vivi_sankey(tech::Tech;time=1,valueIndex=0,heatExergy=false,showHeat=true) = vivi_sankey(tech.i,[tech],tech.o,[],[],time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)
vivi_sankey(problem::Problem;time=1,valueIndex=0,heatExergy=false,showHeat=true) = vivi_sankey(problem.i,problem.p,problem.o,problem.ut,problem.st,time=time,valueIndex=valueIndex,heatExergy=heatExergy,showHeat=showHeat)

# Graphs
function vivi_plot(techs,utils,type;time=1)
    # Preparation
    t_q=[]
    max_segments_per_tech,loads_lims_per_tech = segments_loads_per_tech(techs)
    
    for (τ,tech) in enumerate(techs)        
        size = maximum(tech.f)
        load = tech.f[time]/size
        
        if size !=0 && load !=0
            loads = loads_lims_per_tech[τ]
            index = 1
            while !(loads[index] <= load <=loads[index+1])
                index+=1
            end

            for heat in tech.h
                a,b = get_linear_parameters(loads,heat.pw,max_segments_per_tech[τ])            
                hx = heat.q*(tech.f[time]*a[index]+b[index]*size)
                push!(t_q,Heat(q=hx,Ts=heat.Ts,Tt=heat.Tt))
            end
        end
    end
    techsTq = length(t_q)

    max_segments_per_tech,loads_lims_per_tech = segments_loads_per_tech(utils)
    for (τ,tech) in enumerate(utils)

        size = maximum(tech.f)
        load = tech.f[time]/size
        if size !=0 && load !=0
            loads = loads_lims_per_tech[τ]
            index = 1
            while !(loads[index] <= load <=loads[index+1])
                index+=1
            end

            for heat in tech.h
                a,b = get_linear_parameters(loads,heat.pw,max_segments_per_tech[τ])            
                hx = heat.q*(tech.f[time]*a[index]+b[index]*size)
                push!(t_q,Heat(q=hx,Ts=heat.Ts,Tt=heat.Tt))            
            end
        end
    end

    if type == "icc"
        plot = plot_icc(t_q,techsTq)
    elseif type == "cc"
        plot=plot_cc(t_q)
    elseif type == "gcc"
        plot=plot_gcc(t_q)
    end

    return plot
end

vivi_icc(techs,utils;time=1) = vivi_plot(techs,utils,"icc",time=time)
vivi_icc(techs::Tech;time=1) = vivi_plot([techs],[],"icc",time=time)
vivi_icc(problem::Problem;time=1) = vivi_plot(problem.p,problem.ut,"icc",time=time)

vivi_gcc(techs,utils;time=1) = vivi_plot(techs,utils,"gcc",time=time)
vivi_gcc(techs::Tech;time=1) = vivi_plot([techs],[],"gcc",time=time)
vivi_gcc(problem::Problem;time=1) = vivi_plot(problem.p,problem.ut,"gcc",time=time)

vivi_cc(techs,utils;time=1) = vivi_plot(techs,utils,"cc",time=time)
vivi_cc(techs::Tech;time=1) = vivi_plot([techs],[],"cc",time=time)
vivi_cc(prob::Problem;time=1) = vivi_plot(prob.p,prob.ut,"cc",time=time)