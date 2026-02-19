# Module to plot parameter space 

module ParameterDiagram

using Printf
using LinearAlgebra
using ColorSchemes
using Colors
using UnPack
using GeometryBasics
using Random
using Makie
using CairoMakie
using StaticArrays
using SparseArrays
using CircularArrays
using FromFile
using DrWatson
using GeometryBasics: area

function parameterDiagram(params)

    @unpack γ,
        Λ_00,
        Λ_01,
        Λ_11 = params

    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
    fig2 = Figure(size=(600,600))
    grid = fig2[1,1] = GridLayout()
    ax1 = Axis(grid[1,1],aspect=1)

    ax1.title = "Floppy boundary parameter comparison, Γ=-0.13Λ"
    ax1.xlabel = "Λ"
    ax1.ylabel = "Γ"
    

    Λ_vec = range(-1,0.1,length=200)
    Γ_vec = -0.13 .* Λ_vec
    lines!(ax1,Λ_vec,Γ_vec,color=:red,linewidth=2,label="Floppy boundary: Γ = -0.13 Λ")

    scatter!(ax1,[Λ_00],[γ],color=:blue,markersize=10,label="Λ_00, γ")
    scatter!(ax1,[Λ_01],[γ],color=:green,markersize=10,label="Λ_01, γ")
    scatter!(ax1,[Λ_11],[γ],color=:orange,markersize=10,label="Λ_11, γ")
    axislegend(ax1)

    return fig2 


end # end function 

function P_effsDiagram(t,ax,fig,matrices)

    @unpack P_effs,
        cellAreas = matrices 

    sumP_effA_i = sum(cellAreas.*P_effs)

    scatter!(ax,t,sumP_effA_i,color=:purple,markersize=5)

    return fig

end # end function 

function U_iDiagram(t,ax,fig,matrices,params)

    @unpack cellAreas,
        cellPerimeters,
        edgeLengths,
        Λs = matrices
    @unpack γ = params

    energy1 = sum((1/2).*(cellAreas .- 1).^2)
    energy2 = sum((γ/2).*(cellPerimeters).^2)
    energy3 = sum(Λs .* edgeLengths)
    totalEnergy = energy1+energy2+energy3

    scatter!(ax,t,totalEnergy,color=:purple,markersize=5)

    return fig

end # end funciton

function ξ_iDiagram(t,ax,fig,matrices,params)

    @unpack cellAreas,
        ξs = matrices

    sumAᵢξᵢ = sum(cellAreas .* ξs)

    scatter!(ax,t,sumAᵢξᵢ,color=:purple,markersize=5)

    return fig

end # end function


export parameterDiagram, P_effsDiagram, U_iDiagram, ξ_iDiagram

end # end module 