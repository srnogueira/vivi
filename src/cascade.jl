#=
Heat cascated constraint
=#

"Heat transfer"
mutable struct HeatStruct
    h::Float64
    Ts::Float64
    Tt::Float64
    dtExtra::Float64
    loadEffect::Vector{Vector{Real}}
end

" Array of source shifted tempertures "
function shiftSourceT(stArray,dt)
    shiftHot(st)=st.Ts>st.Tt ? st.Ts-dt-st.dtExtra : st.Ts+st.dtExtra
    return [shiftHot(st) for st in stArray]
end

" Temperature interval array "
function tempInterval(stArray,dt)
    foo=shiftSourceT(stArray,dt)
    unique!(foo)
    sort!(foo,rev=true)
    return foo
end

" Temperature interval array - for plot"
function tempInterval2(stArray,dt)
    bar=zeros(length(stArray)*2)
    i = 1
    for st in stArray
        if st.Ts>st.Tt
            bar[i] = st.Ts-dt/2-st.dtExtra
            bar[i+1] = st.Tt-dt/2-st.dtExtra
        else
            bar[i] = st.Ts+dt/2+st.dtExtra
            bar[i+1] = st.Tt+dt/2+st.dtExtra
        end
        i+=2
    end
    unique!(bar)
    sort!(bar,rev=true)
    return bar
end

" Heat cascade "
function heatCascade(stHX,utHot,utCold;dt=15,forMILP=false,forplot=false,forLP=false)
    # Create the array of streams
    stArray=copy(stHX)
    append!(stArray,utHot)
    append!(stArray,utCold)

    # Temperature intervals
    if forplot
        Tk=tempInterval2(stArray,dt)
    else
        Tk=tempInterval(stArray,dt)
    end

    if forplot
        locHot=falses(length(Tk)-1,length(utHot))
        locCold=falses(length(Tk)-1,length(utCold))
    else
        locHot=zeros(length(Tk)-1,length(utHot))
        locCold=zeros(length(Tk)-1,length(utCold))
    end

    # Count number of hots and colds
    hotN=0
    coldN=0
    for st in stArray
        if st.Ts>st.Tt
            hotN+=1
        else
            coldN+=1
        end
    end

    # Heat per temperature interval
    Qk=zeros(length(Tk)-1)
    Qk_hot=zeros(length(Tk)-1,hotN)
    Qk_cold=zeros(length(Tk)-1,coldN)
    Qk_st=zeros(length(Tk)-1,hotN+coldN)

    hot=0
    cold=0
    for st in stArray # for every stream
        stC=-st.h/(st.Tt-st.Ts) # (+) hot ; (-) cold
        if forplot
            if st.Tt<st.Ts
                hot+=1
                Tt=st.Tt-dt/2-st.dtExtra
                Ts=st.Ts-dt/2-st.dtExtra
            else
                cold+=1
                Tt=st.Tt+dt/2+st.dtExtra
                Ts=st.Ts+dt/2+st.dtExtra
            end
        else
            if st.Tt<st.Ts
                hot+=1
                Tt=st.Tt-dt-st.dtExtra
                Ts=st.Ts-dt-st.dtExtra
            else
                cold+=1
                Tt=st.Tt+st.dtExtra
                Ts=st.Ts+st.dtExtra
            end
        end
        for k=2:length(Tk) # interval: Tk[k-1]-Tk[k]
            up = (Ts>=Tk[k-1]) & (Tt>=Tk[k-1])
            down = (Ts<=Tk[k]) & (Tt<=Tk[k])
            if ~( up | down ) # not out
                if stC>0 # if hot
                    highT = Ts <= Tk[k-1] ? Ts : Tk[k-1]
                    lowT = Tt <= Tk[k] ? Tk[k] : Tt
                else # if cold
                    highT = Tt <= Tk[k-1] ? Tt : Tk[k-1]
                    lowT = Ts <= Tk[k] ? Tk[k] : Ts
                end
                Qk[k-1]+=stC*(highT-lowT)
                if Tt<Ts
                    Qk_hot[k-1,hot]=stC*(highT-lowT)
                else
                    Qk_cold[k-1,cold]=stC*(highT-lowT)
                end
                Qk_st[k-1,hot+cold] = stC*(highT-lowT)
            end
        end
    end

    # Location of hot utilities
    hotId=0
    for st in utHot # for every stream
        hotId+=1
        for k=2:length(Tk) # interval: Tk[k-1]-Tk[k]
            up = (st.Ts-dt-st.dtExtra>=Tk[k-1]) & (st.Tt-dt-st.dtExtra>=Tk[k-1])
            down = (st.Ts-dt-st.dtExtra<=Tk[k]) & (st.Tt-dt-st.dtExtra<=Tk[k])
            if ~( up | down ) # not out
                if (forplot)
                    locHot[k-1,hotId]=true
                    continue
                end
                if st.Tt<st.Ts
                    hot+=1
                    Tt=st.Tt-dt-st.dtExtra
                    Ts=st.Ts-dt-st.dtExtra
                else
                    cold+=1
                    Tt=st.Tt+st.dtExtra
                    Ts=st.Ts+st.dtExtra
                end
                highT = Ts <= Tk[k-1] ? Ts : Tk[k-1]
                lowT = Tt <= Tk[k] ? Tk[k] : Tt
                locHot[k-1,hotId]=highT-lowT #true
            end
        end
    end

    # Location of cold utilities
    coldId=0
    for st in utCold # for every stream
        coldId+=1
        for k=2:length(Tk) # interval: Tk[k-1]-Tk[k]
            up = (st.Ts+st.dtExtra>=Tk[k-1]) & (st.Tt+st.dtExtra>=Tk[k-1])
            down = (st.Ts+st.dtExtra<=Tk[k]) & (st.Tt+st.dtExtra<=Tk[k])
            if ~( up | down ) # not out
                if (forplot)
                    locCold[k-1,coldId]=true
                    continue
                end
                if st.Tt<st.Ts
                    hot+=1
                    Tt=st.Tt-dt-st.dtExtra
                    Ts=st.Ts-dt-st.dtExtra
                else
                    cold+=1
                    Tt=st.Tt+st.dtExtra
                    Ts=st.Ts+st.dtExtra
                end
                highT = Tt <= Tk[k-1] ? Tt : Tk[k-1]
                lowT = Ts <= Tk[k] ? Tk[k] : Ts
                locCold[k-1,coldId]=highT-lowT#true
            end
        end
    end

    # Return info
    if forplot==true
        return Tk,Qk_hot,Qk_cold,locHot,locCold
    elseif forMILP==true
        return Qk,Tk,Qk_hot,Qk_cold,stArray
    elseif forLP==true
        return Qk_st,locHot,locCold
    else
        return Qk,locHot,locCold
    end
end

"Minimal energy requirement"
function mer(stHX,utHot,utCold)
    # mass flow rates and deltahk
    # Prepare
    Qk,locHot,locCold=heatCascade(stHX,utHot,utCold)
    hotN=length(utHot)
    coldN=length(utCold)
    TkN=length(Qk)

    # Define variables
    mer=Model(HiGHS.Optimizer)
    set_silent(mer)
    @variable(mer,hotUTs[i=1:hotN]>=0)
    @variable(mer,coldUTs[i=1:coldN]>=0)
    @variable(mer,R[i=1:TkN]>=0)

    # Define energy balance constraint - vectorized
    @constraint(mer,con[i=1:TkN],R[i]==Qk[i])
    for k=1:TkN
        k > 1 && set_normalized_coefficient(con[k],R[k-1],-1)
        for hot=1:hotN
            if locHot[k,hot] != 0.0
                set_normalized_coefficient(con[k],hotUTs[hot],-1*locHot[k,hot])
            end
        end
        for cold=1:coldN
            if locCold[k,cold] != 0.0
                set_normalized_coefficient(con[k],coldUTs[cold],1*locCold[k,cold])
            end
        end
    end
    @constraint(mer,Rcon,R[TkN]==0)

    # Calculate the exergy costs - Something to avoid crazy results
    # Constrant about the energy transfer from multiple streams

    # Define objective function
    @objective(mer,Min,sum(hotUTs[hot]*(utHot[hot].Ts-utHot[hot].Tt)+coldUTs[cold]*(utCold[cold].Tt-utCold[cold].Ts) for hot=1:hotN for cold=1:coldN))

    # Solve problem
    optimize!(mer)

    # Return optimal values
    merHot=[abs(value(hotUTs[hot])*(utHot[hot].Ts-utHot[hot].Tt)) for hot=1:hotN]
    merCold=[abs(value(coldUTs[cold])*(utCold[cold].Tt-utCold[cold].Ts)) for cold=1:coldN]

    return merHot,merCold
end
