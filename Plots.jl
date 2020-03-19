#=
Plot_Prop:
- Julia version: 1.3.1
- Author: Abd Alazez Ahmed
- Date: 2020-03-18
=#
#- DEFINITIONS -#
ENV["PYTHON"]="/usr/bin/python"
using Plots
pyplot()
range1 = 0
range2 = 30

#- CONSTANSTS -#

DENSITY = 1.225#0.61041
e=0.8
x=range(range1,stop=range2,length=100)

#- FUNCTIONS -#


function calc_K(s,a)
    1/(pi*e*(s^2)/a);
end
function calc_Vtr_min(w,a,k,CD_0)
    sqrt((2/DENSITY)*(w/a)*sqrt(k/CD_0));
end
function calc_Vinf(Vtr,i)
    Vtr-30 .+(3 .*i);
end
function CL_to_Vinf(CL,w,a)
    sqrt(2*w/(CL*DENSITY*a))
end
function Vinf_to_Index(Vinf,Vtr)
    (Vinf-Vtr+30)/3
end
function calc_CL(w,Vinf,a)
    2 .*w./(DENSITY.*(Vinf.^2).*a)
end
function calc_CD(k,CL,CD_0)
    CD_0.+(K.*(CL.^2))
end
function calc_CL_over_CD(CL,CD)
    CL./CD
end
function calc_ThrustReq(CL_over_CD,w)
    w./CL_over_CD
end
function calc_power(ThrustReq,Vinf)
    ThrustReq.*Vinf
end
function calc_RC(Pr,Pa,w)
    (137208.8 .-Pr)./w
end

#- INPUT -$

#$- REPLACE WITH READ FOR CONSOLE USE -$#
print("Weight: ")
w =  parse(Float64,readline())
print("Wing Area: ")
a = parse(Float64,readline())
print("Span: ")
s = parse(Float64,readline())
print("CD_0: ")
CD_0 = parse(Float64,readline())
print("CL_max: ")
CL_MAX = parse(Float64,readline())




#- CALCULATIONS -#

K = calc_K(s,a);
Vmin = calc_Vtr_min(w,a,K,CD_0);
Vinf = calc_Vinf(Vmin,x);
CL = calc_CL(w,Vinf,a);
CD = calc_CD(K,CL,CD_0);
CLoverCD = calc_CL_over_CD(CL,CD);
ThrustReq = calc_ThrustReq(CLoverCD,w);
Power = calc_power(ThrustReq,Vinf);
PowerA = calc_power(16179.9,Vinf);
VinfCrit = CL_to_Vinf(CL_MAX,w,a);
RC = calc_RC(Power,PowerA,w)
index_crit = floor(Integer,Vinf_to_Index(VinfCrit,Vmin)*100/range2-range1)+1+range1;
index_min = floor(Integer,Vinf_to_Index(Vmin,Vmin)*100/range2-range1)+1+range1;



#- PLOTTING FUNCTIONS -#
function plot_ThrustReq_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],ThrustReq[index_crit : length(ThrustReq)],reuse=false,label="Thrust Required",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    hline!([16179.9],line=(:solid,2,:blue))
    plot!(Vinf[1 : index_crit],ThrustReq[1 : index_crit],reuse=false,label="Thrust Required_flase",line=(:dash,3,:red))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("Thrust required")
    xlabel!("V inf")
    return p
end
function plot_Pow_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],Power[index_crit : length(ThrustReq)],reuse=false,label="Power Required",line=(:solid,3,:green))
    plot!(Vinf[index_crit : length(Vinf)],PowerA[index_crit : length(ThrustReq)],reuse=false,label="Power Required",line=(:solid,3,:blue))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    plot!(Vinf[1 : index_crit],Power[1 : index_crit],reuse=false,label="Power Required_flase",line=(:dash,3,:red))
    plot!(Vinf[1 : index_crit],PowerA[1 : index_crit],reuse=false,label="Power Required_flase",line=(:dash,3,:blue))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("Power required")
    xlabel!("V inf")
    return p
end
function plot_RC_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],RC[index_crit : length(RC)],reuse=false,label="RC",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    plot!(Vinf[1 : index_crit],RC[1 : index_crit],reuse=false,label="RC_flase",line=(:dash,3,:red))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("RC")
    xlabel!("V inf")
    return p
end
function plot_ThrustReq_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],ThrustReq[index_crit : length(ThrustReq)],reuse=false,label="Thrust Required",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    hline!([16179.9],line=(:solid,2,:blue))
    plot!(Vinf[1 : index_crit],ThrustReq[1 : index_crit],reuse=false,label="Thrust Required_flase",line=(:dash,3,:red))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("Thrust required")
    xlabel!("V inf")
    return p
end
function plot_Coof_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],CL[index_crit : length(CL)],reuse=false,label="CL",line=(:solid,3,:blue))
    plot!(Vinf[index_crit : length(Vinf)],CD[index_crit : length(CD)],label="CD",line=(:solid,3,:red))
    plot!(Vinf[index_crit : length(Vinf)],CLoverCD[index_crit : length(CLoverCD)],label="CL/CD",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    plot!(Vinf[1 : index_crit],CL[1 : index_crit],reuse=false,label="CL_flase",line=(:dash,3,:blue))
    plot!(Vinf[1 : index_crit],CD[1 : index_crit],label="CD_false",line=(:dash,3,:red))
    plot!(Vinf[1 : index_crit],CLoverCD[1 : index_crit],label="CL/CD_false",line=(:dash,3,:green))
    annotate!(VinfCrit,-1,text(round(VinfCrit,digits=4), 7, :left,:buttom, :red))
    annotate!(Vmin,-1,text(round(Vmin,digits=4), 7, :left, :buttom, :green))
    ylabel!("CL\nCD\nCL/CD")
    xlabel!("V inf")
    return p
end


#- PLOTTING -#

RCPLT = plot_RC_Vinf();
print("Max RC = ")
print(maximum(RC))
print("Vinf Max RC = ")
print(Vinf[findfirst(x-> x == maximum(RC),RC)])

print("ThetaMax = ")
print(atand(maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)])))
angleX=range(0,stop=Vmin,length=100)
plot!(angleX,angleX .* maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)]))
circleX=range(20/sqrt(maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)])^2 +1),stop=20,length=100)
plot!(circleX,sqrt.(400 .- circleX .^ 2))
annotate!(20,0,text((string("Theta max =",atand(maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)])),"degrees")), 7, :left,:buttom, :black))
savefig("RC")
#=
Plot_Jet:
- Julia version: 1.3.1
- Author: Abd Alazez Ahmed
- Date: 2020-03-18
=#
#- DEFINITIONS -#
ENV["PYTHON"]="/usr/bin/python"
using Plots
pyplot()
range1 = 0
range2 = 60

#- CONSTANSTS -#

DENSITY = 1.225#0.61041
e=0.81
x=range(range1,stop=range2,length=100)

#- FUNCTIONS -#


function calc_K(s,a)
    1/(pi*e*(s^2)/a);
end
function calc_Vtr_min(w,a,k,CD_0)
    sqrt((2/DENSITY)*(w/a)*sqrt(k/CD_0));
end
function calc_Vinf(Vtr,i)
    Vtr-30 .+(3 .*i);
end
function CL_to_Vinf(CL,w,a)
    sqrt(2*w/(CL*DENSITY*a))
end
function Vinf_to_Index(Vinf,Vtr)
    (Vinf-Vtr+30)/3
end
function calc_CL(w,Vinf,a)
    2 .*w./(DENSITY.*(Vinf.^2).*a)
end
function calc_CD(k,CL,CD_0)
    CD_0.+(K.*(CL.^2))
end
function calc_CL_over_CD(CL,CD)
    CL./CD
end
function calc_ThrustReq(CL_over_CD,w)
    w./CL_over_CD
end
function calc_power(ThrustReq,Vinf)
    ThrustReq.*Vinf
end
function calc_RC(Pr,Pa,w)
    (Pa .- Pr)./w
end

#- INPUT -$

#$- REPLACE WITH READ FOR CONSOLE USE -$#
print("Weight: ")
w =  parse(Float64,readline())
print("Wing Area: ")
a = parse(Float64,readline())
print("Span: ")
s = parse(Float64,readline())
print("CD_0: ")
CD_0 = parse(Float64,readline())
print("CL_max: ")
CL_MAX = parse(Float64,readline())




#- CALCULATIONS -#

K = calc_K(s,a);
Vmin = calc_Vtr_min(w,a,K,CD_0);
Vinf = calc_Vinf(Vmin,x);
CL = calc_CL(w,Vinf,a);
CD = calc_CD(K,CL,CD_0);
CLoverCD = calc_CL_over_CD(CL,CD);
ThrustReq = calc_ThrustReq(CLoverCD,w);
Power = calc_power(ThrustReq,Vinf);
PowerA = calc_power(32470 ,Vinf);
VinfCrit = CL_to_Vinf(CL_MAX,w,a);
RC = calc_RC(Power,PowerA,w)
index_crit = floor(Integer,Vinf_to_Index(VinfCrit,Vmin)*100/range2-range1)+1+range1;
index_min = floor(Integer,Vinf_to_Index(Vmin,Vmin)*100/range2-range1)+1+range1;



#- PLOTTING FUNCTIONS -#
function plot_ThrustReq_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],ThrustReq[index_crit : length(ThrustReq)],reuse=false,label="Thrust Required",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    hline!([16179.9],line=(:solid,2,:blue))
    plot!(Vinf[1 : index_crit],ThrustReq[1 : index_crit],reuse=false,label="Thrust Required_flase",line=(:dash,3,:red))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("Thrust required")
    xlabel!("V inf")
    return p
end
function plot_Pow_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],Power[index_crit : length(ThrustReq)],reuse=false,label="Power Required",line=(:solid,3,:green))
    plot!(Vinf[index_crit : length(Vinf)],PowerA[index_crit : length(ThrustReq)],reuse=false,label="Power Required",line=(:solid,3,:blue))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    plot!(Vinf[1 : index_crit],Power[1 : index_crit],reuse=false,label="Power Required_flase",line=(:dash,3,:red))
    plot!(Vinf[1 : index_crit],PowerA[1 : index_crit],reuse=false,label="Power Required_flase",line=(:dash,3,:blue))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("Power required")
    xlabel!("V inf")
    return p
end
function plot_RC_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],RC[index_crit : length(RC)],reuse=false,label="RC",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    plot!(Vinf[1 : index_crit],RC[1 : index_crit],reuse=false,label="RC_flase",line=(:dash,3,:red))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("RC")
    xlabel!("V inf")
    return p
end
function plot_ThrustReq_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],ThrustReq[index_crit : length(ThrustReq)],reuse=false,label="Thrust Required",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    hline!([16179.9],line=(:solid,2,:blue))
    plot!(Vinf[1 : index_crit],ThrustReq[1 : index_crit],reuse=false,label="Thrust Required_flase",line=(:dash,3,:red))
    annotate!(VinfCrit,ThrustReq[index_crit]*1.1,text((string("(",round(VinfCrit,digits=4)," , ",round(ThrustReq[index_crit],digits=4),")")), 7, :left,:buttom, :red))
    annotate!(Vmin,ThrustReq[index_min]*1.1,text((string("(",round(Vmin,digits=4)," , ",round(ThrustReq[index_min],digits=4),")")), 7, :left,:buttom, :green))
    ylabel!("Thrust required")
    xlabel!("V inf")
    return p
end
function plot_Coof_Vinf()
    p = plot(Vinf[index_crit : length(Vinf)],CL[index_crit : length(CL)],reuse=false,label="CL",line=(:solid,3,:blue))
    plot!(Vinf[index_crit : length(Vinf)],CD[index_crit : length(CD)],label="CD",line=(:solid,3,:red))
    plot!(Vinf[index_crit : length(Vinf)],CLoverCD[index_crit : length(CLoverCD)],label="CL/CD",line=(:solid,3,:green))
    vline!([VinfCrit],line=(:dot,2,:red))
    vline!([Vmin],line=(:dot,2,:green))
    plot!(Vinf[1 : index_crit],CL[1 : index_crit],reuse=false,label="CL_flase",line=(:dash,3,:blue))
    plot!(Vinf[1 : index_crit],CD[1 : index_crit],label="CD_false",line=(:dash,3,:red))
    plot!(Vinf[1 : index_crit],CLoverCD[1 : index_crit],label="CL/CD_false",line=(:dash,3,:green))
    annotate!(VinfCrit,-1,text(round(VinfCrit,digits=4), 7, :left,:buttom, :red))
    annotate!(Vmin,-1,text(round(Vmin,digits=4), 7, :left, :buttom, :green))
    ylabel!("CL\nCD\nCL/CD")
    xlabel!("V inf")
    return p
end


#- PLOTTING -#

RCPLT = plot_RC_Vinf();
print("Max RC = ")
print(maximum(RC))
print("Vinf Max RC = ")
print(Vinf[findfirst(x-> x == maximum(RC),RC)])
print("ThetaMax = ")
print(atand(maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)])))
angleX=range(0,stop=Vmin,length=100)
plot!(angleX,angleX .* maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)]))
circleX=range(20/sqrt(maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)])^2 +1),stop=20,length=100)
plot!(circleX,sqrt.(400 .- circleX .^ 2))
annotate!(20,0,text((string("Theta max =",atand(maximum(RC[index_crit : length(RC)] ./ Vinf[index_crit : length(Vinf)])),"degrees")), 7, :left,:buttom, :black))
savefig("RC")
    
