using HousingLifeCycleWorkhorse
using BenchmarkTools
using Plots
## Basic usage example. Solves decision problems, then simulates, then plots policy functions, simulated moments, and life-cycle profiles of one individual
par=benchpar(ns=2,nh=4)
par.shgrd[1,:] = [true, true, false, false] # Can only rent smallest two units
par.shgrd[2,:] = [false, true, true, true] # Can't own smallest unit

@time val,pol = solve_decproblems(par);
@time aggmom, pan = simulate(par,pol,val);
#
ih = 1
is = 1
ie = 1
iν = 1
ip = 2
ids = [1] # Includ more id's to show more households

p = plot(
    plot( plot_title="Policy Functions",
        heatmap(par.agegrd[begin:end-10],par.xgrd,pol.c[:,ie,iν,is,ih,ip,begin:end-10],title="c"),
        heatmap(par.agegrd,par.xgrd,pol.b[:,ie,iν,is,ih,ip,:],title="b"),
        heatmap(par.agegrd,par.xgrd,val.V[:,ie,iν,is,ih,ip,1:end-1],title="V"),
        heatmap(par.agegrd,par.xgrd,par.hgrd[pol.ih[:,ie,iν,is,ih,ip,:]].*replace(isequal.(0.0,par.sgrd[pol.is[:,ie,iν,is,ih,ip,:]]),0 .=>-Inf),title="Rent: h",color=cgrad(:viridis,max(length(par.hgrd),2), categorical = true)),
        heatmap(par.agegrd,par.xgrd,par.hgrd[pol.ih[:,ie,iν,is,ih,ip,:]].*replace(isequal.(1.0,par.sgrd[pol.is[:,ie,iν,is,ih,ip,:]]),0 .=>-Inf),title="Own : h",color=cgrad(:viridis,max(length(par.hgrd),2), categorical = true)),
        heatmap(par.agegrd,par.xgrd,par.hgrd[pol.ih[:,ie,iν,is,ih,ip,:]].*replace((0 .< par.sgrd[pol.is[:,ie,iν,is,ih,ip,:]] .< 1),0 .=>-Inf),title="Shared : h",color=cgrad(:viridis,max(length(par.hgrd),2), categorical = true)),
    ),
    plot(plot_title="Simulated Moments",
        plot(aggmom.age,aggmom.x .+ [fill(0,length(par.agegrd)) aggmom.xse -aggmom.xse],title="Wealth",linecolor=[:black :grey70 :grey70]),
        plot(aggmom.age,aggmom.b .+ [fill(0,length(par.agegrd)) aggmom.bse -aggmom.bse],title="Bond",linecolor=[:black :grey70 :grey70]),
        plot(aggmom.age,aggmom.c .+ [fill(0,length(par.agegrd)) aggmom.cse -aggmom.cse],title="Consumption",linecolor=[:black :grey70 :grey70]),
        plot(aggmom.age,aggmom.own*100,title="Own (%)",ylims=(0,100)),
        plot(aggmom.age,aggmom.rent*100,title="Rent (%)",ylims=(0,100)),
        plot(aggmom.age,(1.0 .- (aggmom.own + aggmom.rent))*100,title="Partial (%)",ylims=(0,100)),
        plot(aggmom.age,aggmom.h,title="h",ylims=(par.hgrd[1],par.hgrd[end])),
        plot(aggmom.age,aggmom.ip,title="ip"),
        plot(aggmom.age,aggmom.feasible,title="feasible",ylims=(0.95,1)),
        legend=false
    ),
    plot(plot_title="Simulated Life Cycle Moments",
            scatter(pan.age[1:76],[sh2string.(par.sgrd[pan.is[76*(id-1) + 1:76*id]],par.hgrd[pan.ih[76*(id-1) + 1:76*id]]) for id in ids]
            ,title="House choice"),
            scatter(pan.age[1:76],[pan.x[76*(id-1) + 1:76*id] for id in ids],title="Wealth"),
            scatter(pan.age[1:76],[pan.c[76*(id-1) + 1:76*id] for id in ids],title="Consumption"),
            scatter(pan.age[1:76],[pan.b[76*(id-1) + 1:76*id] for id in ids],title="Bond"),
            scatter(pan.age[1:76],[par.pgrd[pan.ip[76*(id-1) + 1:76*id]] for id in ids],title="Price"),
            scatter(pan.age[1:76],[par.νgrd[pan.iν[76*(id-1) + 1:76*id]] for id in ids],title="Income Shock"),
            legend=false,left_margin = 10Plots.mm
    ),
    layout=grid(3,1),size=(1200,1400).*0.8
)

##
