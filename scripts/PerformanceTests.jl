include("../src/HousingLifeCycleWorkhorse.jl")


## Testing
par=benchpar()
@time val,pol = solve_decproblems(par);
@time aggmom, pan = simulate(par,pol,val);

@btime val,pol = solve_decproblems(par);
@btime aggmom, pan = simulate(par,pol,val);
