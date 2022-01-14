# Julia packages
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using Plots
using UnPack
using DrWatson
@quickactivate

# Local modules
includet("$(projectdir())/src/TopologyChange.jl"); using .TopologyChange
includet("$(projectdir())/src/CreateRunDirectory.jl"); using .CreateRunDirectory
includet("$(projectdir())/src/Visualise.jl"); using .Visualise
includet("$(projectdir())/src/Initialise.jl"); using .Initialise
includet("$(projectdir())/src/Iterate.jl"); using .Iterate
