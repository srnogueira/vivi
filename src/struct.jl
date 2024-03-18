
"""
    Heat(h::Real,Ts::Real,Tt::Real,dT::Real,pw::Matrix{Real})::Heat
    Heat stream representation

    # Arguments:
    - 'h' : Heat trasfer rate [kW]
    - 'Ts' : Source temperature [K]
    - 'Tt' : Target temperature [K]
    - 'dT' : Minimal temperature difference contribution [K]
    - 'pw' : Piece-wise linearization matrix [load_1 factor_1; load_2 factor_2]
    - 'n' : name
"""
Base.@kwdef struct Heat
    q::Real = 0.0
    Ts::Real = 0.0
    Tt::Real = 0.0
    dT::Real = 7.5
    pw::Matrix{Real} = [0 0;1 1]
    n::String = "-"
end

MAXSIZE = 1E3

"""
    ResourceType(n::String,c::Matrix{Real},u::String)
    Represent the characteristics of the resource type

    Arguments:
    - 'n' : Name of the resource
    - 'c' : Cost matrix [time,type]
    - 'u' : Rate unit
"""
Base.@kwdef mutable struct ResourceType
    n::String = "-"
    c::Matrix{Real} = [1;;]
    u::String = "-"
end

"""
    Resource(t::ResourceType,r::Vector{Real},pw::Matrix{Real})
    Represent a resource stream of a certain type

    Arguments:
    - 't' : Type of resource 
    - 'r' : Rate of stream transfer
    - 'pw' : Piece-wise linearization matrix
"""
Base.@kwdef mutable struct Resource
    t::ResourceType
    r::Vector{Real} = [0.0]
    pw::Matrix{Real} = [0 0; 1 1]
end

abstract type ProblemUnit end

"""
    Tech(n::String,i::Vector{Resource},o::Vector{Resource},h::Vector{Heat},c::Vector{Matrix{Real}},s::Tuple{Real,Real},l::Tuple{Real,Real},rp::Tuple{Real,Real},f::Vector{Real},off::Bool)
    Representation of a technology (process or utility)
    
    Arguments:
    - 'n' : Name
    - 'i' : Array of input resources
    - 'o' : Array of output resources
    - 'h' : Array of heat streams
    - 's' : (min, max) sizes
    - 'c' : Vector of costs matrixs [size,cost] for each type. First column is the sizes values
    - 'l' : (min, max) loads
    - 'rp' : (min,max) ramps
    - 'off' : can be switch off? 
    - 'f' : size factor
    - 'f_t' : size factor per time
"""
Base.@kwdef mutable struct Tech <:ProblemUnit
    n::String = "-"
    i::Vector{Resource} = Array{Resource,1}()
    o::Vector{Resource} = Array{Resource,1}()
    h::Vector{Heat} = Array{Heat,1}()
    c::Vector{Matrix{Real}} = [[0 0; MAXSIZE 0]]
    s::Tuple{Real,Real} = (0, MAXSIZE)
    l::Tuple{Real,Real} = (0, 1)
    rp::Tuple{Real,Real} = (-1, 1)
    off::Bool = true
    f::Real = 1
    f_t::Vector{Real} = [1]
end

"""
    Storage(n::String,t::ResourceType,a::Real,c::Vector{Matrix{Real}},s::Tuple{Real,Real},rp::Tuple{Real,Real},f::Real,f_t::Vector{Real})
    Representation of an ideal storage unit

    Arguments:
    - 'n' : Name
    - 't' : Resource type
    - 'a' : Initial amount stored (as a % of strage size)
    - 'c' : Vector of costs matrixs [size,cost] for each type. First column is the sizes values
    - 'rp' : (min,max) ramps
    - 'f' : size factor
    - 'f_t' : size factor per time
"""
Base.@kwdef mutable struct Storage <:ProblemUnit
    n::String = "-"
    t::ResourceType = Resource()
    a::Real = 0
    c::Vector{Matrix{Real}} = [[0 0; MAXSIZE 0]]
    s::Tuple{Real,Real} = (0, MAXSIZE)
    rp::Tuple{Real,Real} = (-1, 1)
    f::Real = 1
    f_t::Vector{Real} = [1]
end

"""
    Problem(i::Vector{Reasource},o::Vector{Resource},p::Vector{Tech},ut::Vector{Tech},st::Vector{Storage})
    Encapsulates information of a problem

    Arguments:
    - 'i' : Input resources
    - 'o' : Output resources
    - 'p' : Tech processes
    - 'ut' : Tech utilities
    - 'st' : Storage utilities
"""
Base.@kwdef mutable struct Problem
    i::Vector{Resource} = Array{Resource,1}()
    o::Vector{Resource} = Array{Resource,1}()
    p::Vector{Tech} = Array{Tech,1}()
    ut::Vector{Tech} = Array{Tech,1}() 
    st::Vector{Storage} = Array{Storage,1}()
end