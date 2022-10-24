# HousingLifeCycleWorkhorse

[![Build Status](https://github.com/eirikbrandsaas/HousingLifeCycleWorkhorse.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eirikbrandsaas/HousingLifeCycleWorkhorse.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/eirikbrandsaas/HousingLifeCycleWorkhorse.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eirikbrandsaas/HousingLifeCycleWorkhorse.jl)

Note: Code is incomplete and preliminary. Missing PDF of model layout.

## A simple life-cycle model of homeownership
1. To use, add the package by URL. Enter the `pkg` mode by typing `]` in an active julia session and then 
`add https://github.com/eirikbrandsaas/HousingLifeCycleWorkhorse.jl` 
    - Or, similarly, `using Pkg;Pkg.add("https://github.com/eirikbrandsaas/HousingLifeCycleWorkhorse.jl")`
2. Run the following lines to solve the model under the benchmark parameters:
 ```julia
##
using HousingLifeCycleWorkhorse

par = benchpar()
V,pol = solve_decproblems(par)
agg,pan = simulate(par,pol,V)
```
3. ...
4. Profit
    
