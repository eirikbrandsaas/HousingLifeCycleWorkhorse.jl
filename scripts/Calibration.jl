include("../src/HousingLifeCycleWorkhorse.jl")

# Issues
# - Start at age 20

## Initialize
par = benchpar(ns=2,nh=4,nx=70,ms=0.04,mb=0.02)
par.shgrd[1,:] = [true, true, false, false] # Can only rent smallest two units
par.shgrd[2,:] = [false, true, true, true] # Can't own smallest unit

@time val,pol = solve_decproblems(par);
@time aggmom, pan = simulate(par,pol,val);

moments = CSV.read(joinpath(@__DIR__,"../data/modelinputs/empirical_moments.csv"),DataFrame)
momg = groupby(moments,:FW)
po = plot(momg[1].age,momg[1].own*100,label="data",ylabel="%")
plot!(po,aggmom.age,aggmom.own*100,label="model")

px = plot(momg[1].age,momg[1].finwealth,label="data",ylabel="1000s 1995 USD")
plot!(px,aggmom.age,aggmom.fx,label="model")

plot(po,px,xlim=(18,50),plot_title="Model Fit",xlabel="Age")
