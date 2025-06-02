
#
#  Stretch.jl
#  VertexModel
#
#  Created by Natasha Cowleu on 02/06/2025.
#
#
# Functions relating to stetching a monolayer

module Stretch

# Julia packages
using FromFile
using LinearAlgebra
using StaticArrays
using UnPack
using FastBroadcast
using SparseArrays
using GeometryBasics


    function stretchCells(R,t, params,matrices)
        @unpack R_initial = matrices
        @unpack λs, stretchType, tStretch = params
        #stretch monolayer, map R_x->(1 + \lambda)R_x, Ry->R-y/(1+\lambda)

        ## This only works for non dividing monolayers currently
        ## interpolate division points?

        #R = reinterpret(SVector{2,Float64}, integrator.u)


        if stretchType=="uniaxial"
            Λ=@SMatrix[
            1+λs 0.0
            0.0 1/(1+λs)
            ]
        elseif stretchType=="biaxial"
            Λ=@SMatrix[
            (1+λs) 0.0
            0.0 (1+λs)
            ]
        end

        stretch=Λ-I(2) #Lambda dot
    
        dRt=([stretch*k for k in R_initial] )./tStretch

        Rt= R_initial+dRt.*t

        R_final=[Λ*k for k in R_initial] 

        return Rt, R_final
    end


    export stretchCells
end