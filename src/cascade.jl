using JuMP, HiGHS
using Plots, DelimitedFiles

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

"""
    c_p(h::Heat)::Real
    Constant pressure heat capacity [kW/K]

    # Arguments
    - 'h' : Heat struct
"""
function c_p(h::Heat)::Real
    return h.q/(h.Ts-h.Tt)
end

"""
    sTs(h::Heat)::Real
    Shifted temperature source [K]

    # Arguments
    - 'h' : Heat struct
"""
function sTs(h::Heat)::Real
    return h.Ts>h.Tt ? h.Ts-h.dT : h.Ts+h.dT
end

"""
    sTt(h::Heat)::Real
    Shifted temperature target [K]

    # Arguments
    - 'h' : Heat struct
"""
function sTt(h::Heat)::Real
    return h.Ts>h.Tt ? h.Tt-h.dT : h.Tt+h.dT
end

"""
    tempInterval(heats::Vector{Heat})::Vector{Real}
    Array of unique and sorted shifted temperatures

    # Arguments
    - 'heats' : Array of heat streams
"""
function tempInterval(heats::Vector{Heat})::Vector{Real}
    temps=[sTs(heat) for heat in heats]
    temps=vcat(temps,[sTt(heat) for heat in heats])
    unique!(temps)
    sort!(temps,rev=true)
    return temps
end

"""
    heatCascade(heats::Vector{Heat})::Matrix{Real}
    Heat table [interval , stream] 

    Arguments:
    - 'heats' : Array of heat streams
"""
function heatCascade(heats::Vector{Heat})::Matrix{Real}
    Tk=tempInterval(heats)

    Qk_st=zeros(length(Tk)-1,length(heats))
    for (i,heat) in enumerate(heats) # for every stream
        for k in eachindex(Tk) # interval: Tk[k-1]-Tk[k]
            k == 1 && continue # skip the first
            up = (sTs(heat)>=Tk[k-1]) & (sTt(heat)>=Tk[k-1])
            down = (sTs(heat)<=Tk[k]) & (sTt(heat)<=Tk[k])
            
            if ~( up | down ) # not out
                if c_p(heat)>0 # if hot
                    highT = sTs(heat) <= Tk[k-1] ? sTs(heat) : Tk[k-1]
                    lowT = sTt(heat) <= Tk[k] ? Tk[k] : sTt(heat)
                else
                    highT = sTt(heat) <= Tk[k-1] ? sTt(heat) : Tk[k-1]
                    lowT = sTs(heat) <= Tk[k] ? Tk[k] : sTs(heat)
                end
                Qk_st[k-1,i] = c_p(heat)*(highT-lowT)
            end
        end
    end

    return Qk_st
end

"""
    separateHotCold(Qk_st::Matrix{Real})::Tuple[Matrix{Real},Matrix{Real}]
    Separates the hot streams from cold streams in the complete heat table
    Hot streams have positive values, while cold streams not.

    Arguments:
    - 'Qk_st' : Matrix of [interval, streams]
"""
function separateHotCold(Qk_st::Matrix{Real})::Tuple{Matrix{Real},Matrix{Real}}
    Qk_hot = Qk_st[:,[sum(Qk_st[:,i])>0 for i in eachindex(Qk_st[1,:])]]
    Qk_cold = Qk_st[:,[sum(Qk_st[:,i])<=0 for i in eachindex(Qk_st[1,:])]]

    return Qk_hot,Qk_cold
end

"""
    mer(streams::Vector{Heat})::Tuple{Real,Real}
    Minimal energy requirement for array of streams in (heating, cooling) format

    Arguments:
    - 'streams' : Array of heat streams
"""
function mer(heats::Vector{Heat})::Tuple{Real,Real}
    Qk = vec(sum(heatCascade(heats),dims=2)) 
    K=length(Qk)

    # Define variables
    mer=Model(HiGHS.Optimizer)
    set_silent(mer)
    @variable(mer,hotUT>=0)
    @variable(mer,coldUT>=0)
    @variable(mer,R[i=1:K]>=0)

    # Define energy balance constraint - vectorized
    @constraint(mer,con[i=1:K],R[i]==Qk[i])
    for k=2:K
        set_normalized_coefficient(con[k],R[k-1],-1)
    end
    set_normalized_coefficient(con[begin],hotUT,-1)
    set_normalized_coefficient(con[end],coldUT,1)

    @constraint(mer,Rcon,R[K]==0)

    # Define objective function
    @objective(mer,Min,hotUT+coldUT)

    # Solve problem
    optimize!(mer)

    # Return optimal values
    merHot=value(hotUT)
    merCold=value(coldUT)

    return (merHot,merCold)
end

"""
    findRange(array::Vector{<:Real})::Tuple{Real,Real}
    Returns the pair of index that are the first unique numbers

    Arguments:
    - 'array': array of values
"""
function findRange(array::Vector{<:Real})::Tuple{Real,Real}
    n=length(array)
    first=1
    for i in 1:n
        if round(array[i],digits=3) !== round(array[i+1],digits=3)
            first=i
            break
        end
    end
    last=n
    for i in n:-1:2
        if round(array[i],digits=3) !== round(array[i-1],digits=3)
            last=i
            break
        end
    end
    return first,last
end


"""
    plot_cc(heats::Vector{Heat};save::Bool = false)
    Plots the hot and cold composite curves

    Arguments:
    - 'heats' : Vector of heats
    - 'save' : Option to save values in a csv file
"""
function plot_cc(heats::Vector{Heat};save::Bool=false)

    Tk = tempInterval(heats)
    Qk_st = heatCascade(heats)
    Qk_hot,Qk_cold = separateHotCold(Qk_st)
    merh,merc = mer(heats)

    cold_cc=vec(sum(Qk_cold,dims=2)) # into vector
    push!(cold_cc,0.0)

    hot_cc=vec(sum(Qk_hot,dims=2))
    push!(hot_cc,0.0)

    # Cascade
    n=length(Tk)-1
    cold_cc = [merc-sum(cold_cc[i:n]) for i=1:n] ./ 1000
    hot_cc = [sum(hot_cc[i:n]) for i=1:n] ./ 1000

    # Removing duplicates (long vertical line)
    ih,oh = findRange(hot_cc)
    ic,oc = findRange(cold_cc)

    # Temperature in celsius
    Tk = Tk .- 273.15

    # Plot
    fig=plot(cold_cc[ic:oc],Tk[ic:oc],label="Cold",lw=2)
    plot!(fig,hot_cc[ih:oh],Tk[ih:oh],label="Hot",lw=2)

    if save
        writedlm("hot.csv",[hot_cc[ih:oh] Tk[ih:oh]],",")
        writedlm("cold.csv",[cold_cc[ic:oc] Tk[ic:oc]],",")
    end

    xlabel!("Net heat load [MW]")
    ylabel!("\$T_{shifted}\$ [°C]")
    return fig
end

"""
    compute_gcc(heats::Vector{Heat};utility::Bool=false)::Tuple{Vector{Real},Vector{Real}}

    Delivers the x,y points for the gcc plot

    Arguments:
    - 'heats' : Array of heat streams
    - 'utility' : If the heat streams should be treated as utilities for icc plots
"""
function compute_gcc(heats::Vector{Heat};utility::Bool=false)::Tuple{Vector{Real},Vector{Real}}
    # Get values
    Tk = tempInterval(heats)
    Qk_st = heatCascade(heats)
    Qk_hot,Qk_cold = separateHotCold(Qk_st)
    merh,merc = mer(heats)
    cold_cc=vec(sum(Qk_cold,dims=2)) # into vector
    hot_cc=vec(sum(Qk_hot,dims=2)) # into vector

    # Get composite curve
    R=zeros(length(Tk))
    signal = utility ? -1 : 1
    R[1]=(sum(merh))/1000
    for i in 2:length(Tk)
        R[i]=R[i-1]+signal*(hot_cc[i-1]+cold_cc[i-1])/1000
    end

    # Get the interval
    x,xo = findRange(R)

    # Temperature in celsius
    Tk = Tk .- 273.15

    return R[x:xo],Tk[x:xo]
end

"""
    plot_gcc(heats::Vector{Heat};save::Bool=false)
    Plots the grand composite curves

    Arguments:
    - 'heats' : Array of heat streams
    - 'save' : Option to save values in a csv file
"""
function plot_gcc(heats::Vector{Heat};save::Bool=false)
    R,Tk = compute_gcc(heats)

    # Plot GCC
    save && writedlm(filename,[R[x:xo] Tk[x:xo]],",")
    fig=plot(R,Tk,label="",lw=2,xlims=(0,Inf))
    xlabel!("Net heat load [MW]")
    ylabel!("\$T_{shifted}\$ [°C]")

    return fig
end

"""
    plot_icc(heats::Vector{Heat},processes::Int64;save::Bool=false)
    Plots the integrated composite curves
    * Assumes that heats is ordered as processes streams before utility streams

    Arguments:
    - 'heats' : Array of heat streams
    - 'processes' : Number of streams attributed to processes
    - 'save' : Option to save values in a csv file
"""
function plot_icc(heats::Vector{Heat},processes::Int64;save::Bool=false)
    # Get values    
    R_p,Tk_p = compute_gcc(heats[1:processes])
    R_u,Tk_u = compute_gcc(heats[(processes+1):end],utility=true)

    # Correcting R_u
    R_u = R_u .- R_u[1] .+ R_p[1]

    fig=plot(R_p,Tk_p,label="Processes",lw=2)
    fig=plot!(R_u,Tk_u,label="Utilities",lw=2,xlims=(min(R_u...)-0.01,Inf))

    xlabel!("Net heat load [MW]")
    ylabel!("\$T_{shifted}\$ [°C]")

    if save
        writedlm("icc_process.csv",[R_p Tk_p],",")
        writedlm("icc_utility.csv",[R_u Tk_u],",")
    end

    return fig
end

"""
    test_cascade()

    Simple test function for the minimal energy requirement problem
"""
function test_cascade()
    # Test problem 4SP1 from Papoulias and Grossmann (1983) 
    C1 = Heat(q=762,Ts=60+273,Tt=160+273,dT=5)
    H2 = Heat(q=589,Ts=160+273,Tt=93+273,dT=5)
    C3 = Heat(q=876,Ts=116+273,Tt=260+273,dT=5)
    H4 = Heat(q=1171,Ts=249+273,Tt=138+273,dT=5)

    heats = [C1,H2,C3,H4]
    merHot,merCold = mer(heats)

    println("Computed solution is $(round(merHot)) kW and $(round(merCold)) kW for hot and cold utilities")
    println("Reported solution is 128 kW and 250 kW for hot and cold utilities")

    fig_cc = plot_cc(heats)
    fig_gcc = plot_gcc(heats)

    n_p = length(heats)
    HUT = Heat(q=merHot,Ts=270+273,Tt=269+273,dT=5)
    CUT = Heat(q=merCold,Ts=38+273,Tt=82+273,dT=5)
    
    heats = vcat(heats,[HUT,CUT])
    fig_icc = plot_icc(heats,n_p)

    return plot(fig_cc,fig_gcc,fig_icc,layout=(1,3))
end