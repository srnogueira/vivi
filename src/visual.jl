using Plots, DelimitedFiles
theme(:dao) # Pretty plots

include("graphs.jl")

function plot_cc(hx_array,h,c,merh,merc)
    # Get heat exchanged by intervals
    Tk,Qk_hot,Qk_cold,locHot,locCold= heatCascade(hx_array,h,c,forplot=true)
    cold_cc=vec(sum(Qk_cold,dims=2)) # into vector
    ΔH_cold=sum(merc)
    ΔH_hot=sum(merh)
    push!(cold_cc,0.0)

    hot_cc=vec(sum(Qk_hot,dims=2))
    push!(hot_cc,0.0)

    # Cascade
    n=length(Tk)-1
    cold_cc=-cold_cc
    for i in 1:n
        cold_cc[i]=(sum(cold_cc[i:n])+ΔH_cold)/1000
    end

    for i in 1:n
        hot_cc[i]=(sum(hot_cc[i:n]))/1000
    end

    # Removing duplicates (long vertical line)
    ih=1
    for i in 1:n
        if round(hot_cc[i],digits=3) !== round(hot_cc[i+1],digits=3)
            ih=i
            break
        end
    end

    ic=1
    for i in 1:n
        if round(cold_cc[i],digits=3) !== round(cold_cc[i+1],digits=3)
            ic=i
            break
        end
    end

    oh=n
    for i in n:-1:2
        if round(hot_cc[i],digits=3) !== round(hot_cc[i-1],digits=3)
            oh=i
            break
        end
    end

    oc=n
    for i in n:-1:2
        if round(cold_cc[i],digits=3) !== round(cold_cc[i-1],digits=3)
            oc=i
            break
        end
    end

    # Temperature in celsius
    for i = 1:length(Tk)
        Tk[i] -= 273.15
    end

    # Plot
    fig=plot(cold_cc[ic:oc],Tk[ic:oc],label="Cold CC",lw=2,legend=:topleft,linecolor=:cadetblue)
    plot!(fig,hot_cc[ih:oh],Tk[ih:oh],label="Hot CC",lw=2,legend=:topleft,linecolor=:indianred)

    #writedlm("hot.txt",[hot_cc[ih:oh] Tk[ih:oh]],"\t")
    #writedlm("cold.txt",[cold_cc[ic:oc] Tk[ic:oc]],"\t")

    xlabel!("Heat load [MW]")
    ylabel!("\$T_{shifted}\$ [°C]")
    return fig
end

# Does not have multiple streams
function plot_gcc(hx_array,h,c,merh,merc)
    # Get values
    Tk,Qk_hot,Qk_cold,locHot,locCold= heatCascade(hx_array,h,c,forplot=true)
    cold_cc=vec(sum(Qk_cold,dims=2)) # into vector
    hot_cc=vec(sum(Qk_hot,dims=2)) # into vector

    # Get composite curve
    R=zeros(length(Tk))
    R[1]=(sum(merh))/1000 #hot_cc[1]+cold_cc[1]+ #DOES NOT WORK IF HOTUT is lower than process hots
    for i in 2:length(Tk)
        R[i]=R[i-1]+(hot_cc[i-1]+cold_cc[i-1])/1000
    end

    # Get the interval
    x=1
    for i in 1:length(Tk)
        if R[i] !== R[i+1]
            x=i
            break
        end
    end
    n=length(Tk)

    xo=1
    for i in length(Tk):-1:2
        if R[i] !== R[i-1]
            xo=i
            break
        end
    end

    # Temperature in celsius
    for i = 1:length(Tk)
        Tk[i] -= 273.15
    end

    # Plot GCC
    fig=plot(R[x:xo],Tk[x:xo],label="GCC",lw=2,xlims=(0,Inf))

    #=
    # Utilities - one line
    if merc[1] !== 0.0
        Ts = c[1].Ts-273.15
        Tt = c[1].Tt-273.15
        plot!(fig,[0,merc[1]/1000],[Ts,Tt],label="Cold UT",lw=3,xlims=(0,Inf))
    end

    # don't know how to determine if the curve should be after or before the previous one
    #plot!(fig,[merc[1]/1000,(merc[1]+merc[3])/1000],[c[3].Ts,c[3].Tt],label="Cold UT2",lw=3,xlims=(200,205),ylims=(200,400))

    if merh[1] !== 0.0
        Ts = h[1].Ts-273.15
        Tt = h[1].Tt-273.15
        plot!(fig,[0,merh[1]/1000],[Ts,Tt],label="Hot UT",lw=3,xlims=(0,Inf))
    end
    =#
    xlabel!("Heat load [MW]")
    ylabel!("\$T_{shifted}\$ [°C]")

    #writedlm("gcc.txt",[R[x:xo] Tk[x:xo]],"\t")
    return fig
end

function plot_icc(hx_array,merh,merc,processes)
    # Get values
    h = [HeatStruct(0.0,3000.0,2999.0)]
    c = [HeatStruct(0.0,1.0,2.0)]
    Tk,Qk_hot,Qk_cold,locHot,locCold= heatCascade(hx_array[1:processes],h,c,forplot=true)
    cold_cc=vec(sum(Qk_cold,dims=2)) # into vector
    hot_cc=vec(sum(Qk_hot,dims=2)) # into vector

    # Get composite curve
    R=zeros(length(Tk))
    R[1]=(sum(merh))/1000 #hot_cc[1]+cold_cc[1]+ #DOES NOT WORK IF HOTUT is lower than process hots
    for i in 2:length(Tk)
        R[i]=R[i-1]+(hot_cc[i-1]+cold_cc[i-1])/1000
    end

    # Get the interval
    x=1
    for i in 1:length(Tk)
        if R[i] !== R[i+1]
            x=i
            break
        end
    end
    n=length(Tk)

    xo=1
    for i in length(Tk):-1:2
        if R[i] !== R[i-1]
            xo=i
            break
        end
    end

    # Temperature in celsius
    for i = 1:length(Tk)
        Tk[i] -= 273.15
    end

    # Plot GCC
    fig=plot(R[x:xo],Tk[x:xo],label="Processes",lw=2,xlims=(0,Inf))

    # Utilities
    Tk,Qk_hot,Qk_cold,locHot,locCold= heatCascade(hx_array[(processes+1):end],h,c,forplot=true)
    cold_cc=vec(sum(Qk_cold,dims=2)) # into vector
    hot_cc=vec(sum(Qk_hot,dims=2)) # into vector

    # Get composite curve
    R=zeros(length(Tk))
    R[1]=(sum(merh))/1000 #hot_cc[1]+cold_cc[1]+ #DOES NOT WORK IF HOTUT is lower than process hots
    for i in 2:length(Tk)
        R[i]=R[i-1]-(hot_cc[i-1]+cold_cc[i-1])/1000
    end

    # Get the interval
    x=1
    for i in 1:length(Tk)
        if R[i] !== R[i+1]
            x=i
            break
        end
    end
    n=length(Tk)

    xo=1
    for i in length(Tk):-1:2
        if R[i] !== R[i-1]
            xo=i
            break
        end
    end

    # Temperature in celsius
    for i = 1:length(Tk)
        Tk[i] -= 273.15
    end

    # Plot GCC
    fig=plot!(R[x:xo],Tk[x:xo],label="Utilities",lw=2,xlims=(min(R...)-0.01,Inf))

    xlabel!("Heat load [MW]")
    ylabel!("\$T_{shifted}\$ [°C]")

    #writedlm("icc.txt",[R[x:xo] Tk[x:xo]],"\t")
    return fig
end


function vivi_plot(techs,utils,type)
    # Preparation
    t_q=[]
    for tech in techs
        aux = deepcopy(tech.heat) # required
        for heat in aux
            heat.h*=tech.size[1]
        end
        append!(t_q,aux)
    end
    techsTq = length(t_q)
    for tech in utils
        aux = deepcopy(tech.heat) # required
        for heat in aux
            heat.h*=tech.size[1]
        end
        append!(t_q,aux)
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

vivi_icc(techs,utils) = vivi_plot(techs,utils,"icc")
vivi_gcc(techs,utils) = vivi_plot(techs,utils,"gcc")
vivi_cc(techs,utils) = vivi_plot(techs,utils,"cc")
