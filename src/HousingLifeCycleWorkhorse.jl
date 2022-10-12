module HousingLifeCycleWorkhorse

## Load packages
using BenchmarkTools
using Interpolations
using Optim
using DelimitedFiles
using StableRNGs
using Statistics
using Distributions
using PanelDataTools # Easy lead/lag/spell functions
using DataFrames
using Plots

## Load modules/functions
include("convenience.jl") # Functions not part of the model (numerical helper functions)
include("fundamentals.jl") # Model fundamentals (budget equation, law of motion, utility...)
include("decisionproblem.jl") # All code for solving HH dec. problem
include("simulation.jl") # Simulate the households
include("sharedownership.jl") # Functions that are specific to shared ownership (i.e., s ∈(0,1) instead of s ∈ {0,1})

## Export of main program things
export benchpar
export solve_decproblems
export simulate

## export of model fundamentals
export util, LoM, r

## Export of various smaller utilities
export sh2string
export rescale_parameters

end
