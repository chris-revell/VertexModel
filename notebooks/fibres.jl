using CairoMakie
using StaticArrays
using Random
using Distributions
using LinearAlgebra
using GeometryBasics

nFibres = 100
lFibre = 300.0

xRange = (0.0,1000.0)
yRange = (0.0,1000.0)



for i=1:nFibres
    pos = rand(Uniform(xRange...),2)
    orientation = normalize!(rand(Uniform(-1.0,1.0),2))
    pts = Point2([pos.-lFibre.*orientation./2, pos.+lFibre.*orientation./2])
    lines!(ax,)