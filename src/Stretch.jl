
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

@from "VertexModelContainers.jl" using VertexModelContainers


    function stretchFunc(t, params)
        #stretch monolayer, map R(0)_x->(1 + \lambda)R(0)_x, R(0)_y->R(0)_y/(1+\lambda)

        ## This only works for non dividing monolayers currently
        ## interpolate division points?
        @unpack λs, stretchType, tStretch = params

        if stretchType=="uniaxial"
            Λt=Diagonal([1+(λs*t)/tStretch, 1-(λs*t)/(tStretch*(1+λs))])
        elseif stretchType=="biaxial"
            Λt=Diagonal([1+(λs*t)/tStretch, 1+(λs*t)/tStretch])
        end
        return Λt
    end

    function dstretchFunc(t, params)

        @unpack λs, stretchType, tStretch = params

        if stretchType=="uniaxial"
            dΛt=Diagonal([λs/tStretch, -λs/(tStretch*(1+λs))])
        elseif stretchType=="biaxial"
           dΛt=Diagonal([λs/tStretch, λs/tStretch])
        end

        return dΛt
    end

    function stretchCells(R,t, params,matrices)
        @unpack R_membrane, Rt, R_final= matrices
        @unpack λs, stretchType, tStretch, tMemChange, nVerts = params
        #stretch monolayer, map R_x->(1 + \lambda)R_x, Ry->R-y/(1+\lambda)

        ## This only works for non dividing monolayers currently
        ## interpolate division points?



        # if stretchType=="uniaxial"
        #     Λt=Diagonal([1+(λs*t)/tStretch, 1/(1-(λs*t)/(tStretch*(1+λs)))])
        # elseif stretchType=="biaxial"
        #     Λt=Diagonal([1+(λs*t)/tStretch, 1+(λs*t)/tStretch])
        # end
        # stretch=Λ-I(2) #Lambda dot

        #dRt=([stretch*k for k in R_membrane] )./tStretch
        
        # Rt= R_membrane+dRt.*t

        #R_final=[Λ*k for k in R_membrane] 

        #@show tMemChange

        # ΛtMemt=(inv(stretchFunc(tMemChange, params))*stretchFunc(t, params)) #Λ(t)
        # #ΛtMemτ=(inv(stretchFunc(tMemChange, params))*stretchFunc(tStretch, params)) #Λ(tStretch)

        # dΛtMemt=(inv(stretchFunc(tMemChange, params))*dstretchFunc(t, params)) #dΛ(t)
       
        

        # dRt=[dΛtMemt*kk for kk in R_membrane]
        # Rt.=[ΛtMemt*kk for kk in R_membrane]

        # if t<=tStretch
        #     R_final.=[ΛtMemτ*kk for kk in R_membrane]
        # else
        #     R_final.=R_membrane
        # end


        #Try having new vertices not stuck down, so stretch force only acting on original points

        #So don't update R_membrane
        #Do stretch on original verts eg 1:nVerts, remaining dRt ->0, remaining Rt, R_final -> R
        n=size(R_membrane)[1]

        dRt=fill(SVector{2,Float64}(zeros(2)), nVerts)
            

        dRt[1:n].=[dstretchFunc(t, params)*kk for kk in R_membrane]
        Rt[1:n].=[stretchFunc(t, params)*kk for kk in R_membrane]
        R_final[1:n].=[stretchFunc(tStretch, params)*kk for kk in R_membrane]

        Rt[n+1:end].=R[n+1:end]
        R_final[n+1:end].=R[n+1:end]

        
        return dRt
    end


    export stretchFunc, dstretchFunc,stretchCells
end