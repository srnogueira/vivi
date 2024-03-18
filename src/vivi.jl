using JuMP

include("struct.jl")
include("cascade.jl")   # Pinch functions
include("graphs.jl")

"""
    add_tech_sizes!(model::Model,techs::Vector{<:ProblemUnit})
    adds the size of techs and min/max size constraints

    Arguments:
    - 'model' : JuMP model
    - 'techs' : Array of Techs in the problem
"""
function add_tech_sizes!(model::Model,techs::Vector{<:ProblemUnit})::Nothing
    n_techs = length(techs)
    @variable(model,f[τ=1:n_techs]) 
    @constraint(model,[τ=1:n_techs],f[τ]>=techs[τ].s[1])
    @constraint(model,[τ=1:n_techs],f[τ]<=techs[τ].s[2])
    return nothing
end

"""
    linear_coeffs(points::Matrix{<:Real})
    returns the a, b coefficients from a matrix of points

    Arguments:
    - 'points': matrix of piecewise linear points
"""
function linear_coeffs(points::Matrix{<:Real})::Tuple{Vector{Real},Vector{Real}}
    n_seg = length(points[:,1])-1
    a = zeros(n_seg)
    b = zeros(n_seg)
    for s=1:n_seg
        a[s] = (points[s+1,2]-points[s,2])/(points[s+1,1]-points[s,1])
        b[s] = points[s+1,2]-a[s]*points[s+1,1]
    end
    return a,b
end

"""
    create_capex_MILP!(model::Model,techs::Vector{<:ProblemUnit};valueIndex::Int64=1)
    adds piecewise linearization of CAPEX

    Arguments:
    - 'model' : JuMP model
    - 'techs' : Array of Techs in the problem
    - 'valueIndex' : Index of cost vector
"""
function create_capex_MILP!(model::Model,techs::Vector{<:ProblemUnit};valueIndex::Int64=1)::Nothing
    # Check number of segments for each Tech
    n_techs = length(techs)
    n_segments = Vector{Int}(undef,n_techs)

    for τ in eachindex(techs)
        points = techs[τ].c[valueIndex]
        n_segments[τ] = length(points[:,1])-1 # Two points for every segment
        if n_segments[τ] <= 0
            throw("Tech $(techs[τ].n) does not have cost points")
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

    # Add segment sizes
    @variable(model,f_s[τ=1:n_techs,s=1:maximum(n_segments); matrix[τ,s] > 0]>=0)
    @constraint(model,[τ=1:n_techs],sum(f_s[τ,:])==f[τ])
    for τ in eachindex(techs)
        points = techs[τ].c[valueIndex]
        for s=1:n_segments[τ]
            @constraint(model,points[s,1]*y_s[τ,s]<=f_s[τ,s])
            @constraint(model,f_s[τ,s]<=points[s+1,1]*y_s[τ,s])
        end
    end

    # Create capex expression
    a = zeros(n_techs,maximum(n_segments))
    b = zeros(n_techs,maximum(n_segments))
    for τ in eachindex(techs)
        points = techs[τ].c[valueIndex]
        a[τ,:],b[τ,:] = linear_coeffs(points)
    end

    @expression(model,capex,sum(f_s[τ,s]*a[τ,s]+y_s[τ,s]*b[τ,s] for (τ,s) in eachindex(y_s)))
    return nothing
end

"""
    create_capex_LP!(model::Model,techs::Vector{<:ProblemUnit};valueIndex::Int64=1)::nothing
    adds linear approximation of CAPEX

    Arguments:
    - 'model' : JuMP model
    - 'techs' : Array of Techs in the problem
    - 'valueIndex' : Index of cost vector
"""
function create_capex_LP!(model::Model,techs::Vector{<:ProblemUnit};valueIndex::Int64=1)::Nothing
    n_techs = length(techs)

    # Create continuous variable
    add_tech_sizes!(model,techs)
    f = model[:f]

    @expression(model,capex,sum(f[τ]*techs[τ].c[valueIndex][end,2]/techs[τ].c[valueIndex][end,1] for τ= 1:n_techs))
    return nothing
end

"""
    load_limits(techs::Vector{Tech})::Vector{Real}
    returns info about the segments load limits from all resources and heats in a tech

    Arguments:
    - 'techs' : Array of Techs in the problem
"""
function load_limits(tech::Tech)::Vector{Real}
    loads = Set()
    for i in vcat(tech.i,tech.o,tech.h)
        for point in i.pw[:,1]
            tech.l[1] <= point <= tech.l[2] && push!(loads,point)
        end
    end

    if tech.l[1] == tech.l[2]
        loads_lims = [tech.l[1],tech.l[2]]
    else
        push!(loads,tech.l[1])
        push!(loads,tech.l[2])
        loads_lims = sort(collect(loads))
    end
    return  loads_lims
end

"""
    get_linear_parameters(loads::Vector{<:Real},points::Matrix{Float64},max_segments::Int64)
    returns the slope and intersection parameters of each linear segment for Resource and Heat

    Arguments:
    - 'tech': Tech 
    - 'points': matrix of piecewise linearization of resource/heat
"""
function get_linear_parameters(tech::Tech,points::Matrix{Real})::Tuple{Vector{Real},Vector{Real}}
    loads = load_limits(tech)
    max_segments = length(loads)-1
    a = zeros(max_segments)
    b = zeros(max_segments)
    if !isempty(points)
        # Calculate segment coefficients from provided points
        a_points,b_points = linear_coeffs(points)
        # Match coefficients with the tech segments
        p = 1
        for s = 1:max_segments
            if loads[s] == points[p,1] # Found a match
                a[s] = a_points[p]
                b[s] = b_points[p]
                p+=1
            else
                if (points[p,1] <= loads[s] <= points[p+1,1])
                    # (A) is contained in a larger segment
                    a[s] =  a_points[p] 
                    b[s] = b_points[p]
                else
                    # (B) is undefined
                    throw("Undefined load in piecewise linearization of $(tech.n)")
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
    create_tech_reformulation!(model::Model,techs::Vector{Tech},resources::Vector{Resources},n_time::Int64)
    creates matrixes with the tech resources and heat expressions using the reformulation strategy

    Arguments:
    - 'model' : a JuMP model
    - 'techs' : Array of techs
    - 'resources' Array of resources
    - 'n_time' : number of time points 
"""
function create_tech_reformulation!(model::Model,techs::Vector{Tech},resources::Vector{String},n_time::Int64)
    # Find number of segments for each resource and heat of each tech
    n_techs = length(techs)
    loads_lims_per_tech = [load_limits(t) for t in techs]
    max_segments_per_tech = [length(l)-1 for l in loads_lims_per_tech]
    max_segments = maximum(max_segments_per_tech)

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
        tech.off == false && @constraint(model,[t=1:n_time],sum(δ[τ,:,t])==1) # Can't be shut down
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
    r_ts .= @expression(model,0)

    for (τ,tech) in enumerate(techs)
        n_in = length(tech.i)
        for (i,res) in enumerate(vcat(tech.i,tech.o))
            signal = i > n_in ? 1 : -1
            r = findfirst(x -> x == res.t.n,resources)
            # Finding piewise representation
            a,b = get_linear_parameters(tech,res.pw)
            # Add resource balance
            if length(res.r) == 1 # Maybe it does not make so much sense
                r_ts[r,:] = [@expression(model,r_ts[r,t]+signal*res.r[1]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ])) for t=1:n_time]
            else
                r_ts[r,:] = [@expression(model,r_ts[r,t]+signal*res.r[t]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ])) for t=1:n_time]
            end
        end
    end

    # Heats
    heats = Array{Heat,1}()
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

    q_kt = Matrix{Any}(undef,n_Tk,n_time)
    q_kt .= @expression(model,0)
 
    start = 0
    for (τ,tech) in enumerate(techs)
        for (h,heat) in enumerate(tech.h)
            Qki = Qk[:,h+start]
            a,b = get_linear_parameters(tech,heat.pw)
            q_kt = [@expression(model,q_kt[k,t]-Qki[k]*sum(ξ[τ,s,t]*b[s]+f_st[τ,s,t]*a[s] for s=1:max_segments_per_tech[τ])) for k in 1:n_Tk, t=1:n_time]
        end
        start += length(tech.h)
    end
    return r_ts,q_kt
end

"""
    create_tech_LP!(model::Model,techs::Vector{Tech},resources::Vector{Resource},n_time::Int64)
    creates matrixes with the tech resources and heat expressions using the simplified strategy

    Arguments:
    - 'model' : a JuMP model
    - 'techs' : Array of techs
    - 'resources' Array of resources
    - 'n_time' : number of time points 
"""
function create_tech_LP!(model::Model,techs::Vector{Tech},resources::Vector{String},n_time::Int64)
    # Find number of segments for each resource and heat of each tech
    n_techs = length(techs)

    # Create the variables - Reformulation strategy (f_ts are defined outside)
    f = model[:f]
    f_t = model[:f_t]
    
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] <= techs[τ].s[2])
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] >= techs[τ].s[1])
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] <= f[τ]*techs[τ].l[2])
    @constraint(model,[t=1:n_time,τ=1:n_techs],f_t[τ,t] >= f[τ]*techs[τ].l[1])

    # Mass contributions
    n_resources = length(resources)
    
    r_ts = Matrix{Any}(undef,n_resources,n_time)
    r_ts .= @expression(model,0)

    for (τ,tech) in enumerate(techs)
        n_in = length(tech.i)
        for (i,res) in enumerate(vcat(tech.i,tech.o))
            signal = i > n_in ? 1 : -1
            r = findfirst(x -> x == res.t.n,resources)
            # Add resource balance
            if length(res.r) == 1 # Maybe it does not make so much sense
                r_ts[r,:] = [@expression(model,r_ts[r,t]+signal*res.r[1]*f_t[τ,t]) for t=1:n_time]
            else
                r_ts[r,:] = [@expression(model,r_ts[r,t]+signal*res.r[t]*f_t[τ,t]) for t=1:n_time]
            end
        end
    end

    # Heats
    heats = Array{Heat,1}()
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
 
    start = 0
    for (τ,tech) in enumerate(techs)
        for h in eachindex(tech.h)
            Qki = Qk[:,h+start]
            q_kt = [@expression(model,q_kt[k,t]-Qki[k]*f_t[τ,t]) for k in 1:n_Tk, t=1:n_time]
        end
        start += length(tech.h)
    end

    return r_ts,q_kt
end

"""
    add_ramping!(model::Model,techs::Vector{Tech},n_time::Int64)
    add ramping constraints for the techs in the optimization model

    Arguments:
    - 'model' : a JuMP model
    - 'techs' : Array of techs
    - 'resources' Array of resources
    - 'n_time' : number of time points 
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

    Argument:
    -'problem' : Problem structure
"""
function number_of_time_points(problem::Problem)
    n_time = 0
    for resource in vcat(problem.i,problem.o)
        n_time = max(n_time,length(resource.r))
    end
    for tech in vcat(problem.p,problem.ut)
        for r in vcat(tech.i,tech.o)
            n_time = max(n_time,length(r.r))
        end
    end
    return n_time
end


"""
    list_of_resources(problem::Problem)
    returns the resources in a problem as a Set

    Arguments:
    - 'problem' : Problem structure
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

"""
    create_model(problem::Problem,solver::DataType;valueIndex::Int64=1,capex=false,LP=false)::Model
    Creates the optimization problem

    Arguments:
    - 'problem' : problem structure
    - 'solver' : a JuMP solver
    - 'valueIndex' : index of value
    - 'capex' : include capex in objective function
    - 'LP' : linear formulation of the problem
"""
function create_model(problem::Problem,solver::DataType;valueIndex::Int64=1,capex::Bool=false,LP::Bool=false)::Model
    # Get info #################################################################
    model = Model(solver)

    # Define variables
    n_time = number_of_time_points(problem)
    resources = list_of_resources(problem)
    techs = vcat(problem.p,problem.ut)
    n_res = length(resources)
    n_techs = length(techs)
    n_store = length(problem.st)

    # Technology size/load variables ##############################################
    if capex
        if LP
            create_capex_LP!(model,vcat(techs,problem.st),valueIndex=valueIndex)
        else
            create_capex_MILP!(model,vcat(techs,problem.st),valueIndex=valueIndex)
        end
    else
        add_tech_sizes!(model,vcat(techs,problem.st))
    end
    f = model[:f]
    @variable(model,f_t[τ=1:n_techs+n_store,t=1:n_time]>=0) # Store

    # Reformulation strategy
    if LP
        r_ts,q_kt = create_tech_LP!(model,techs,resources,n_time) # heat expressions are also created here
    else
        r_ts,q_kt = create_tech_reformulation!(model,techs,resources,n_time) # heat expressions are also created here
    end

    # Resources constraints #######################################################
    @constraint(model,mass[i=1:n_res,t=1:n_time],r_ts[i,t]==0)

    binary_inOut = zeros(n_res,2)
    for input in problem.i
        r = findfirst(x -> x == input.t.n,resources)
        binary_inOut[r,1] = 1
    end
    for output in problem.o
        r = findfirst(x -> x == output.t.n,resources)
        binary_inOut[r,2] = 1
    end
    
    @variable(model,resIN[r=1:n_res,1:n_time;binary_inOut[r,1]==1]>=0)
    @variable(model,resOUT[r=1:n_res,1:n_time;binary_inOut[r,2]==1]>=0)
    
    for r=1:n_res # OBS: got problmes when trying to make it shorter
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

    valIN = zeros(n_res,n_time)
    valOUT = zeros(n_res,n_time)

    for input in problem.i
        i = findfirst(x -> x == input.t.n,resources)
        if length(input.t.c[:,valueIndex]) == n_time
            valIN[i,:] = input.t.c[:,valueIndex]
        elseif length(input.t.c[:,valueIndex]) == 1
            valIN[i,:] = [input.t.c[1,valueIndex] for i =1:n_time]
        else
            throw("mismatched lenghts of arrays")
        end
        
        if length(input.r) == n_time
            for t = 1:n_time
                if input.r[t] != Inf
                    @constraint(model,resIN[i,t] == input.r[t])
                end
            end
        elseif length(input.r) == 1
            if input.r[1] != Inf
                for t = 1:n_time
                    @constraint(model,resIN[i,t] == input.r[1])
                end
            end
        else
            throw("mismatched lenghts of arrays")
        end
    end

    # Outputs
    for output in problem.o
        i = findfirst(x -> x == output.t.n,resources)
        if length(output.t.c[:,valueIndex]) == n_time
            valOUT[i,:] = output.t.c[:,valueIndex]
        elseif length(output.t.c[:,valueIndex]) == 1
            valOUT[i,:] = [output.t.c[1,valueIndex] for i =1:n_time]
        else
            throw("mismatched lenghts of arrays")
        end

        if length(output.r) == n_time
            for t = 1:n_time
                if output.r[t] != Inf
                    @constraint(model,resOUT[i,t] == output.r[t])
                end
            end
        elseif length(output.r) == 1
            if output.r[1] != Inf
                for t = 1:n_time
                    @constraint(model,resOUT[i,t] == output.r[1])
                end
            end
        else
            throw("mismatched lenghts of arrays")
        end
    end
    
    # Storage #########################################################
    # Intial storage constraint
    if n_store >= 1
        @constraint(model,[i=1:n_store,t=1:n_time],f_t[n_techs+i,t]<=problem.st[i].s[2])
        @constraint(model,[i=1:n_store,t=1:n_time],f_t[n_techs+i,t]>=problem.st[i].s[1])
        @constraint(model,[i=1:n_store],f_t[n_techs+i,n_time] == f[n_techs+i]*problem.st[i].a)
        @constraint(model,size2[i=1:n_store,t=1:n_time],f[n_techs+i]>=f_t[n_techs+i,t]) # It should be higher or equal to the highest load
    end

    # Balance of storage
    for (s_index,s) in enumerate(problem.st)
        i = findfirst(x -> x == s.t.n,resources)
        set_normalized_coefficient(mass[i,1],f[n_techs+s_index],s.a)
        set_normalized_coefficient(mass[i,1],f_t[n_techs+s_index,1],-1)
        for t = 2:n_time
            set_normalized_coefficient(mass[i,t],f_t[n_techs+s_index,t-1],1)
            set_normalized_coefficient(mass[i,t],f_t[n_techs+s_index,t],-1)
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
        @objective(model,Max,cost_out-cost_in)
    end

    return model
end

"""
    create_problem_answer(problem::Problem,model::Model)::Problem
    Returns the answer as a Problem

    Arguments:
    - 'problem' : problem structure
    - 'model' : a JuMP model
"""
function create_problem_answer(problem::Problem,model::Model)::Problem
    f_t = model[:f_t]
    f = model[:f]
    resIN = model[:resIN]
    resOUT = model[:resOUT]
    n_techs = length(vcat(problem.p,problem.ut))
    n_time = number_of_time_points(problem)
    n_store = length(problem.st)
    resources = list_of_resources(problem)
    
    inputs_ans = deepcopy(problem.i)
    outputs_ans = deepcopy(problem.o)
    processes_ans = deepcopy(problem.p)
    utils_ans = deepcopy(problem.ut)
    store_ans = deepcopy(problem.st)

    gamma_opt = [round(value(f_t[τ,t]),sigdigits=5) for τ=1:n_techs for t=1:n_time]
    gamma_opt = reshape(gamma_opt,(n_time,n_techs))

    store_opt = [round(value(f_t[s,t]),sigdigits=5) for s=n_techs+1:n_techs+n_store for t=1:n_time]
    store_opt = reshape(store_opt,(n_time,n_store))
    
    # Inputs
    for input in inputs_ans
        input.r=zeros(n_time)
        r_index = findfirst(x -> x == input.t.n,resources)
        for t = 1:n_time
            ans_in = value(resIN[r_index,t])
            input.r[t] = ans_in 
        end
    end

    # Output
    for output in outputs_ans
        output.r=zeros(n_time)        
        r_index = findfirst(x -> x == output.t.n,resources)
        for t = 1:n_time
            ans_out = value(resOUT[r_index,t])
            output.r[t] = ans_out
        end
    end

    # Techs
    for (i,ans) in enumerate(processes_ans)
        ans.f = value(f[i])
        ans.f_t = gamma_opt[:,i]
    end

    # Utilities
    for (i,ans) in enumerate(utils_ans)
        ans.f = value(f[i+length(problem.p)])
        ans.f_t = gamma_opt[:,i+length(problem.p)]
    end

    # Storage
    for (i,ans) in enumerate(store_ans)
        ans.f = value(f[i+n_techs])
        ans.f_t = store_opt[:,i]
    end

    return Problem(i=inputs_ans,p=processes_ans,o=outputs_ans,ut=utils_ans,st=store_ans)
end

"""
    printSolution(model;valueIndex::Int64=1)::Nothing
    Prints the main solution info in the terminal

    Arguments:
    - 'model' : the JuMP model
"""
function printSolution(model)::Nothing
    f = model[:f]
    
    println(solution_summary(model))

    larger = length("Tech")
    larger_2 = length("Size")
    for (i,tech) in enumerate(vcat(problem.p,problem.ut,problem.st))
        larger = length(tech.n) > larger ? length(tech.n) : larger
        larger_2 = length(string(value(f[i]))) > larger_2 ? length(string(value(f[i]))) : larger_2
    end
    
    head = " Tech"
    line = "====="
    subline = "-----"
    while length(head) < larger + 3
        head = head * " "
        line = line * "="
        subline = subline * "-"
    end
    while length(line) < larger + 9 + larger_2
        line = line * "="
        subline = subline * "-"
    end

    println("$head | Size factor")
    println(line)

    for (i,tech) in enumerate(problem.p)
        n = " " * tech.n
        while length(n) < larger + 3
            n = n * " "
        end
        v = value(f[i])
        println("$n | $v")
    end

    println(subline)
    n_p = length(problem.p)
    for (i,tech) in enumerate(problem.ut)
        n = " " * tech.n
        while length(n) < larger + 3
            n = n * " "
        end
        v = value(f[i+n_p])
        println("$n | $v")
    end

    println(subline)
    n_techs = length(vcat(problem.p,problem.ut))
    for (s,store) in enumerate(problem.st)
        n = " " * store.n
        while length(n) < larger + 3
            n = n * " "
        end
        v = value(f[s+n_techs])
        println("$n | $v")
    end
    return nothing
end

"""
    vivi(problem::Problem,solver::DataType;valueIndex::Int64=1,capex=false,LP=false)::Model
    A wrapper to create and solve the problem

    Arguments:
    - 'problem' : problem structure
    - 'solver' : a JuMP solver
    - 'valueIndex' : index of value
    - 'capex' : include capex in objective function
    - 'LP' : linear formulation of the problem
"""
function vivi(problem,solver;valueIndex::Int64=1,capex::Bool=false,LP::Bool=false,print=true)::Problem
    model = create_model(problem,solver,valueIndex=valueIndex,capex=capex,LP=LP)
    optimize!(model)
    print && printSolution(model)
    return create_problem_answer(problem,model)
end