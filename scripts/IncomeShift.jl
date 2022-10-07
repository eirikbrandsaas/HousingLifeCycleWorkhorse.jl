include("../src/HousingLifeCycleWorkhorse.jl")

# Issues
# - look at how much you have to increase high FW to match data

## Initialize
moments = CSV.read(joinpath(@__DIR__,"../data/modelinputs/empirical_moments.csv"),DataFrame)
momg = groupby(moments,:FW)

par = benchpar(ns=2,nh=3,nx=50,ms=0.04,mb=0.02,β=0.97,χ=.066,ne=3,np=3,d=0.1,nν=7,nϵ=2)
par.shgrd[1,:] = [true, false, false] # Can only rent smallest two units
par.shgrd[2,:] = [false, true, true] # Can't own smallest unit

@time val,pol = solve_decproblems(par);
@time aggmom, pan = simulate(par,pol,val);

#
po = plot(ylim=(0,100))
plot!(po,momg[1].age,momg[1].own*100,label="Data (low FW)",ylabel="%",color=:blue,title="Homeownership")
plot!(po,aggmom.age,aggmom.own*100,label="Model",color=:blue,linestyle=:dash)

px = plot(ylim=(0,150))
plot!(px,momg[1].age,momg[1].finwealth,label="data",ylabel="1000s 1995 USD",color=:blue,title="Financial Wealth")
plot!(px,aggmom.age,aggmom.fx,label="model",color=:blue,linestyle=:dash,legend=false)

pb = plot(po,px,xlim=(18,42),plot_title="Model Fit (Low FW)",xlabel="Age",legend=:topleft)
savefig(pb,"tabfig/MatchingMoments.pdf")
pb
##
highFWinc = 1.22
par = benchpar(ns=length(par.sgrd),nh=length(par.hgrd),nx=length(par.xgrd),ms=par.ms,mb=par.mb,β=par.β,χ=par.χ,ne=length(par.egrd),np=length(par.egrd),d=par.d,nν=length(par.νgrd),nϵ=par.nϵ,
    highFWinc=highFWinc)

par.shgrd[1,:] = [true, false, false] # Can only rent smallest two units
par.shgrd[2,:] = [false, true, true] # Can't own smallest unit

@time val,pol = solve_decproblems(par);
@time aggmom, pan = simulate(par,pol,val);

plot!(po,momg[1].age,momg[2].own*100,label="Data (high FW)",ylabel="%",color=:red)
plot!(po,aggmom.age,aggmom.own*100,label="Model (low FW, shifted income)",color=:red,linestyle=:dash)

plot!(px,momg[1].age,momg[2].finwealth,label="data",ylabel="1000s 1995 USD",color=:red)
plot!(px,aggmom.age,aggmom.fx,label="model",color=:red,linestyle=:dot)

plot(po,px,xlim=(18,42),plot_title="Model Fit",xlabel="Age")

p = plot(po,xlim=(18,42),plot_title="Income shift of $(highFWinc)",xlabel="Age",legend=:topleft)
savefig(p,"tabfig/IncomeShift.pdf")
p
