#
#  InitialHexagons.jl
#  VertexModel
#
#  Created by Natasha Cowley on 05/01/2024.
#
#
# Function to create an initial large random system

module RandomInitialSystem

# Julia packages
using LinearAlgebra
using SparseArrays
using StaticArrays
using DrWatson
using FromFile
using DelaunayTriangulation
using StaticArrays
using FromFile
using DelimitedFiles
using Statistics

@from "SenseCheck.jl" using SenseCheck

function cellInitialSystem()
    matdir="C:\\Users\\v35431nc\\Documents\\Lab_Stuff\\Movies_to_track\\100cells\\20151125_1_GSV_GFPtub-CheHis_uu_0p5_MP_fr361\\2024-02-03_15-00-15\\Matrices"
    Rfile=datadir(matdir,"20151125_1_fr361_Matrix_R.txt")
    Afile=datadir(matdir,"20151125_1_fr361_Matrix_A.txt")
    Bfile=datadir(matdir,"20151125_1_fr361_Matrix_B.txt")
    Rcells= readdlm(Rfile, ' ', Float64)
    Acells= readdlm(Afile, ' ', Float64)
    Bcells= readdlm(Bfile, ' ', Float64)

    B=dropzeros!(sparse(Bcells))
    A=dropzeros!(sparse(Acells))
    R=[SVector(Rcells[k,1]/25, Rcells[k,2]/25) for k in 1:size(Rcells)[1]]
    R.-=mean(R, dims=1)

    senseCheck(A, B; marker="read cells")

    
    return A, B, R

end

export cellInitialSystem 

end