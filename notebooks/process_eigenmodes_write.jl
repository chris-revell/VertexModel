s# Import Julia packages
using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using StaticArrays
using BlockArrays
using CairoMakie
using UnPack
using FromFile
using GeometryBasics
using Random
using Colors
using JLD2
using LaTeXStrings
using Glob
using Printf
using ColorSchemes


@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/CellProperties.jl" using CellProperties

folder="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims/new_energy/Hex_relax/relaxed"
files=Glob.glob("new_energy/Hex_relax/relaxed/systemData*Gamma_0.5*.jld2","C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims")
mkpath(datadir(folder,"cell_plots"))
mkpath(datadir(folder,"spectra"))
mkpath(datadir(folder,"cell_modes"))
mkpath(datadir(folder,"vertex_modes"))

mkpath(datadir(folder,"data"))
for f in files
    print(f)
    @unpack R, matrices, params = load(f)
    @unpack nCells,nVerts, γ, L₀ = params
    @unpack cellAreas,cellPerimeters, cellTensions, cellPressures, cellEdgeCount, edgeTangents,vertexAreas,A,Ā,B,B̄,C, ϵ = matrices


    print(' ')
    print(nCells)
    print('\n')

    G=makeGnew(params, matrices) #update for new model
    M=makeM(matrices)
    Lv=makeEvLvnew(M, G, vertexAreas)
    Lc=makeEvLcnew(M, G, vertexAreas)
    X=makeX(params, matrices)
    g=vcat(cellPressures, cellTensions)
    
    gX=Matrix{SMatrix{2,2,Float64,4}}(undef,nVerts,nVerts)
    fill!(gX,@SMatrix zeros(2,2))
    
    for α=1:2*nCells
         gX+=g[α]X[α, :,:]
    end

    dAdr=-1/2*B*Diagonal([(ϵ*T) for T in edgeTangents])*Ā
    dLdr= B̄*Diagonal((edgeTangents)./norm.(edgeTangents))*A

    Mtemp=vcat(dAdr,dLdr)
    Mflat=transpose(reshape(reinterpret(Float64,[vcat(Matrix(Mtemp)[x,:]...) for x in collect(1:2*nCells)]), ( 2*nVerts, 2*nCells)))
        
    E=Diagonal(vertexAreas)
    Eflat=kron(E,Matrix(1.0I, 2, 2))
    
    LcE=Matrix(Lc)
    evalLc,evecLctemp=LAPACK.syev!('V','U',deepcopy(sqrt(G)*Matrix(LcE)*inv(sqrt(G))))
    evecLc=inv(sqrt(G))*evecLctemp
    
    LvE=Matrix(Lv)
    evalLv,evecLvtemp=LAPACK.syev!('V','U',deepcopy(sqrt(Eflat)*Matrix(mortar(LvE))*inv(sqrt(Eflat))))
    evecLv=inv(sqrt(Eflat))*evecLvtemp
    
    H=Matrix(mortar(LvE)).+Matrix(mortar(inv(E)*gX))
    evalH,evecHtemp=LAPACK.syev!('V','U',deepcopy(sqrt(Eflat)*H*inv(sqrt(Eflat))))
    evecH=inv(sqrt(Eflat))*evecHtemp
    
    #Ntemp=Matrix(sqrt.(G)*M*inv(sqrt.(E)))
    # Ntempflat=transpose(reshape(reinterpret(Float64,[vcat(Matrix(Ntemp)[x,:]...) for x in collect(1:2*nCells)]), ( 2*nVerts, 2*nCells)))
    Nflat=sqrt.(G)*Mflat*inv(sqrt(Eflat))
    UNF,sNF,VNF=svd(Nflat, full=true)
    
    Y=inv(sqrt.(G))*UNF
    Z=inv(sqrt.(Eflat))*VNF
    
    D=evecLv'*Matrix(mortar(gX))*evecLv
    DD=D + Diagonal(evalLv)
    evalDD,evecDD=LAPACK.syev!('V','U',deepcopy(DD))

    evecmap=[evecDD[:,x]'*(evecLv'*Matrix(mortar(gX))*evecLv+Diagonal(evalLv))*evecDD[:,x] for x in 1:2nVerts]
    evmapLv=[evecDD[:,x]'*(Diagonal(evalLv))*evecDD[:,x] for x in 1:2nVerts]
    evmapgX=[evecDD[:,x]'*(evecLv'*Matrix(mortar(gX))*evecLv)*evecDD[:,x] for x in 1:2nVerts]



    dat_dir=mkpath(datadir(folder,"data", "Γ_"*string(params.γ)*"_L0_"*string(params.L₀)))

    writedlm(datadir(dat_dir,"LcE_eigenvalues.csv"), evalLc, ',') 
    writedlm(datadir(dat_dir,"LcE_eigenvectors.csv"), evecLc, ',') 
    writedlm(datadir(dat_dir,"LvE_eigenvalues.csv"), evalLv, ',') 
    writedlm(datadir(dat_dir,"LvE_eigenvectors.csv"), evecLv, ',') 
    writedlm(datadir(dat_dir,"HE_eigenvalues.csv"), evalH, ',') 
    writedlm(datadir(dat_dir,"HE_eigenvectors.csv"), evecH, ',') 
    writedlm(datadir(dat_dir,"DppE_tot_diag_sigma_eigenvalues.csv"), evalDD, ',') 
    writedlm(datadir(dat_dir,"DppE_tot_diag_sigma_eigenvectors.csv"), evecDD, ',') 
    writedlm(datadir(dat_dir,"svd_N_Y.csv"), Y, ',') 
    writedlm(datadir(dat_dir,"svd_N_Z.csv"), Z, ',') 
    writedlm(datadir(dat_dir,"svd_N_sigma.csv"), sNF, ',') 
    writedlm(datadir(dat_dir,"evmap.csv"), evecmap, ',') 
    writedlm(datadir(dat_dir,"evmapLv.csv"), evmapLv, ',') 
    writedlm(datadir(dat_dir,"evmapgX.csv"), evmapgX, ',') 

   


end