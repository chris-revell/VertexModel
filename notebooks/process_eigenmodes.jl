# Import Julia packages
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


@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/CellProperties.jl" using CellProperties

function makeE(matrices, params)


    @unpack edgeTangents,A, Ā, B, C= matrices
    @unpack nCells, nEdges, nVerts = params



    vertexAreas=zeros(nVerts)
    edgeMidpointLinks=fill(SVector{2, Float64}(zeros(2)), (nCells, nVerts))

    nzC = findnz(C)
    ikPairs = tuple.(nzC[1],nzC[2])
    for (i,k) in ikPairs
        for j=1:nEdges
            edgeMidpointLinks[i,k] = edgeMidpointLinks[i,k] .+ 0.5.*B[i,j]*edgeTangents[j]*Ā[j,k]
        end
    end


    ##Check k_is length ==1

    for k=1:nVerts

        k_is = findall(x->x!=0, C[:,k])
        if length(k_is) == 1
            edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
        elseif length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
            edgesSharedBy_i2_And_k = findall(x->x!=0, B[k_is[2],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] += 0.5^3*norm([edgeTangents[edgesSharedBy_i2_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i2_And_k[2]]...,0.0])
        else
            vertexAreas[k] = 0.5*norm([edgeMidpointLinks[k_is[1], k]...,0.0]×[edgeMidpointLinks[k_is[2],k]...,0.0])
        end
    end

    E=Diagonal(vertexAreas)

    return(E)

end

folders=Glob.glob("relax_100_cells_fluid/*L₀*","C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims" )

for f in folders
    print(f)
    @unpack R, matrices, params = load(datadir(f,"frameData/systemDataFinal.jld2"))
    @unpack B, Bᵀ, C, cellPositions, cellAreas, cellPerimeters, cellPressures, cellTensions = matrices
    @unpack nCells,nVerts = params

    mkpath(datadir(f,"eigenmodes"))
    print(' ')
    print(nCells)
    print('\n')

    E=makeE(matrices, params)
    Eflat=kron(E,Matrix(1.0I, 2, 2))
    G=makeG(params)
    M=makeM(matrices)
    Lv=makeEvLv(M, G)
    X=makeX(params, matrices)
    g=vcat(cellPressures, -cellTensions)

    gX=Matrix{SMatrix{2,2,Float64,4}}(undef,nVerts,nVerts)
    fill!(gX,@SMatrix zeros(2,2))

    for α=1:2*nCells
        gX+=g[α]X[α, :,:]
    end
    
    LcE=Matrix(tr.(M*(inv(E)*M'))*G)
    evalLc,evecLctemp=LAPACK.syev!('V','U',deepcopy(sqrt(G)*Matrix(LcE)*inv(sqrt(G))))
    evecLc=inv(sqrt(G))*evecLctemp
    
    LvE=inv(E)*Lv
    evalLv,evecLvtemp=LAPACK.syev!('V','U',deepcopy(sqrt(Eflat)*Matrix(mortar(LvE))*inv(sqrt(Eflat))))
    evecLv=inv(sqrt(Eflat))*evecLvtemp
    
    H=Matrix(mortar(LvE)).+Matrix(mortar(inv(E)*gX))
    evalH,evecHtemp=LAPACK.syev!('V','U',deepcopy(sqrt(Eflat)*H*inv(sqrt(Eflat))))
    evecH=inv(sqrt(Eflat))*evecHtemp
    
    Ntemp=Matrix(sqrt.(G)*M*inv(sqrt.(E)))
    Ntempflat=transpose(reshape(reinterpret(Float64,[vcat(Matrix(Ntemp)[x,:]...) for x in collect(1:2*nCells)]), ( 2*nVerts, 2*nCells)))
    UNF,sNF,VNF=svd(Ntempflat, full=true)
    
    Y=inv(sqrt.(G))*UNF
    Z=inv(sqrt.(Eflat))*VNF
    
    D=evecLv'*Matrix(mortar(gX))*evecLv
    DD=D + Diagonal(evalLv)
    evalDD,evecDD=LAPACK.syev!('V','U',deepcopy(DD))

    evecmap=[evecDD[:,x]'*Eflat*(evecLv'*Matrix(mortar(gX))*evecLv+Diagonal(evalLv))*evecDD[:,x] for x in 1:2nVerts]
    evmapLv=[evecDD[:,x]'*Eflat*(Diagonal(evalLv))*evecDD[:,x] for x in 1:2nVerts]
    evmapgX=[evecDD[:,x]'*Eflat*(evecLv'*Matrix(mortar(gX))*evecLv)*evecDD[:,x] for x in 1:2nVerts]

    writedlm(datadir(f,"eigenmodes","LcE_eigenvalues.csv"),evalLc, ',') 
    writedlm(datadir(f,"eigenmodes","LcE_eigenvectors.csv"),evecLc, ',') 
    writedlm(datadir(f,"eigenmodes","LvE_eigenvalues.csv"),evalLv, ',') 
    writedlm(datadir(f,"eigenmodes","LvE_eigenvectors.csv"), evecLv, ',') 
    writedlm(datadir(f,"eigenmodes","HE_eigenvalues.csv"), evalH, ',') 
    writedlm(datadir(f,"eigenmodes","HE_eigenvectors.csv"), evecH, ',') 
    writedlm(datadir(f,"eigenmodes","DppE_tot_diag_sigma_eigenvalues.csv"), evalDD, ',') 
    writedlm(datadir(f,"eigenmodes","DppE_tot_diag_sigma_eigenvectors.csv"), evecDD, ',') 
    writedlm(datadir(f,"eigenmodes","svd_N_Y.csv"), Y, ',') 
    writedlm(datadir(f,"eigenmodes","svd_N_Z.csv"), Z, ',') 
    writedlm(datadir(f,"eigenmodes","svd_N_sigma.csv"), sNF, ',') 

    nv=LinRange(1, 2*nVerts, 2*nVerts)
    nc=LinRange(2*nVerts-2*nCells+1, 2*nVerts, 2*nCells)

    fig = Figure()
    ax=Axis(fig[1, 1], xlabel="n", ylabel="log₁₀λₙ", title="Γ = "*string(params.γ)*", L₀ = "*string(params.L₀))
    
    scatter!(ax,nc, log10.(sort(abs.(sNF.^2))), color=:black,markersize=4, label="σ²ₙ, N")
    scatter!(ax,nv, log10.(sort(abs.(eigvals(H)))), color=:blue,markersize=4, label="λₙ, E⁻¹MᵀGM+E⁻¹gX")
    scatter!(ax,nv, log10.(sort(abs.(evalDD))), color=:red,markersize=4, label="λₙ D+diag(σ²ₙ)")#
    
    #vlines!(ax,2*nVerts+1-((2*nCells)-1), color=:red)
    fig[1, 2] = Legend(fig, ax, framevisible = false)
    save(datadir(f,"eigenmodes","compare_eigenvalues_E_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)

    cellPolygons = makeCellPolygons(R,params,matrices)

    for n=1:2*nCells

        Aevlims=(-maximum(abs.(Y[1:nCells,n])), maximum(abs.(Y[1:nCells, n])))
        Levlims=(-maximum(abs.(Y[ nCells+1:2*nCells,n])), maximum(abs.(Y[nCells+1:2*nCells, n])))
        set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica", fontsize=19)
        fig = Figure(resolution=(1500,500))
    
        a1=Axis(fig[1,1],aspect=DataAspect())
        a2=Axis(fig[2,1],aspect=DataAspect())
        hidedecorations!(a1)
        hidespines!(a1)
        hidedecorations!(a2)
        hidespines!(a2)
        for i=1:nCells
            poly!(a1,cellPolygons[i],color=[Y[1:nCells,n][i]],colormap=:bwr,colorrange=Aevlims, strokecolor=(:black,1.0),strokewidth=1)
            poly!(a2,cellPolygons[i],color=[Y[nCells+1:2*nCells,n][i]],colormap=:bwr,colorrange=Levlims, strokecolor=(:black,1.0),strokewidth=1)
        end
        Label(fig[2,1,Bottom()],"λ_"*string(2*nCells-(n-1))*" = "*@sprintf("%.5E", sNF[n]^2),fontsize = 32)
    
        #hidedecorations!(ax22)
        #hidespines!(ax22)
    
        colsize!(fig.layout,1,Aspect(1,1.0))
    
    
        Colorbar(fig[1,2],limits=colorrange=Aevlims,colormap=:bwr,flipaxis=true)
        Colorbar(fig[2,2],limits=colorrange=Levlims,colormap=:bwr,flipaxis=true)
    
    
        Label(fig[1,1,Left()],string(L"Area"),fontsize = 32, rotation=π/2)
        Label(fig[2,1,Left()],string(L"Perimeter"),fontsize = 32, rotation=π/2)
        #Label( fig[0,:],"Γ = "*string(params.γ)*", L₀ = "*string(params.L₀)*", δL = "*string(params.δL),fontsize = 32, color = (:black, 1))
        resize_to_layout!(fig)
        save(datadir(f,"eigenmodes","eigenmodes_Y$(@sprintf("%03d", 2*nCells-(n-1))).png"),fig)
    
        #display(fig)
    end

    Aevlims=(minimum(abs.(cellAreas[1:nCells])), maximum(abs.(cellAreas[1:nCells])))
    Levlims=(minimum(abs.(cellPerimeters[1:nCells])), maximum(abs.(cellPerimeters[1:nCells])))
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica", fontsize=19)
    fig = Figure(resolution=(1500,500))

    a1=Axis(fig[1,1],aspect=DataAspect())
    a2=Axis(fig[2,1],aspect=DataAspect())
    hidedecorations!(a1)
    hidespines!(a1)
    hidedecorations!(a2)
    hidespines!(a2)
    for i=1:nCells
        poly!(a1,cellPolygons[i],color=cellAreas[i],colormap=:viridis,colorrange=Aevlims, strokecolor=(:black,1.0),strokewidth=1)
        poly!(a2,cellPolygons[i],color=cellPerimeters[i],colormap=:viridis,colorrange=Levlims, strokecolor=(:black,1.0),strokewidth=1)
    end
    #Label(fig[2,1,Bottom()],"λ_"*string(n)*" = "*@sprintf("%.5E", evals[n]),fontsize = 32)

    #hidedecorations!(ax22)
    #hidespines!(ax22)

    colsize!(fig.layout,1,Aspect(1,1.0))


    Colorbar(fig[1,2],limits=colorrange=Aevlims,colormap=:viridis,flipaxis=true)
    Colorbar(fig[2,2],limits=colorrange=Levlims,colormap=:viridis,flipaxis=true)


    Label(fig[1,1,Bottom()],string(L"Area"),fontsize = 32, rotation=0)
    Label(fig[2,1,Bottom()],string(L"Perimeter"),fontsize = 32, rotation=0)
    #Label( fig[0,:],"Γ = "*string(params.γ)*", L₀ = "*string(params.L₀)*", δL = "*string(params.δL),fontsize = 32, color = (:black, 1))
    resize_to_layout!(fig)
    save(datadir(f,"eigenmodes","Area_Perimeter.png"),fig)

    cellShapeTensors=makeShapeTensors(R,params, matrices)
    circularity=getCircularity(params, cellShapeTensors)
    shapeParameter=cellPerimeters./sqrt.(cellAreas)

    Aevlims=(minimum(shapeParameter[1:nCells]), maximum(shapeParameter[1:nCells]))
    Levlims=(minimum(circularity[1:nCells]), maximum(circularity[1:nCells]))
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica", fontsize=19)
    fig = Figure(resolution=(1500,500))

    a1=Axis(fig[1,1],aspect=DataAspect())
    a2=Axis(fig[2,1],aspect=DataAspect())
    hidedecorations!(a1)
    hidespines!(a1)
    hidedecorations!(a2)
    hidespines!(a2)
    for i=1:nCells
        poly!(a1,cellPolygons[i],color=circularity[i],colormap=:viridis,colorrange=Levlims, strokecolor=(:black,1.0),strokewidth=1)
        poly!(a2,cellPolygons[i],color=shapeParameter[i],colormap=:cividis,colorrange=Aevlims, strokecolor=(:black,1.0),strokewidth=1)
    end
    #Label(fig[2,1,Bottom()],"λ_"*string(n)*" = "*@sprintf("%.5E", evals[n]),fontsize = 32)

    #hidedecorations!(ax22)
    #hidespines!(ax22)

    colsize!(fig.layout,1,Aspect(1,1.0))


    Colorbar(fig[1,2],limits=colorrange=Levlims,colormap=:viridis,flipaxis=true)
    Colorbar(fig[2,2],limits=colorrange=Aevlims,colormap=:cividis,flipaxis=true)


    Label(fig[1,1,Bottom()],"Circularity",fontsize = 26, rotation=0)
    Label(fig[2,1,Bottom()],"Shape parameter",fontsize = 26, rotation=0)
    #Label( fig[0,:],"Γ = "*string(params.γ)*", L₀ = "*string(params.L₀)*", δL = "*string(params.δL),fontsize = 32, color = (:black, 1))
    resize_to_layout!(fig)
    save(datadir(f,"eigenmodes","circularity_shape.png"),fig)

    Peff=getPeff(params, matrices)
    cellQ, cellJ=makeCellQandJ(params, matrices)
    cellShearStress=getShearStress(params, matrices, cellJ)

    Plims=(minimum(cellPressures[1:nCells]), maximum(cellPressures[1:nCells]))
    Tlims=(minimum(-cellTensions[1:nCells]), maximum(-cellTensions[1:nCells]))
    Pefflims=(-maximum(abs.(Peff)), maximum(abs.(Peff)))
    ShStlims=(minimum(cellShearStress[1:nCells]), maximum(cellShearStress[1:nCells]))
    set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica", fontsize=19)
    fig = Figure(resolution=(1500,500))

    a11=Axis(fig[1,1],aspect=DataAspect())
    a21=Axis(fig[2,1],aspect=DataAspect())
    a13=Axis(fig[1,3],aspect=DataAspect())
    a23=Axis(fig[2,3],aspect=DataAspect())
    hidedecorations!(a11)
    hidespines!(a11)
    hidedecorations!(a21)
    hidespines!(a21)
    hidedecorations!(a13)
    hidespines!(a13)
    hidedecorations!(a23)
    hidespines!(a23)
    for i=1:nCells
        poly!(a11,cellPolygons[i],color=cellPressures[i],colormap=cgrad(:Blues_9, rev=true),colorrange=Plims, strokecolor=(:black,1.0),strokewidth=1)
        poly!(a21,cellPolygons[i],color=-cellTensions[i],colormap=:plasma,colorrange=Tlims, strokecolor=(:black,1.0),strokewidth=1)
        poly!(a13,cellPolygons[i],color=Peff[i],colormap=:bwr,colorrange=Pefflims, strokecolor=(:black,1.0),strokewidth=1)
        poly!(a23,cellPolygons[i],color=cellShearStress[i],colormap=:plasma,colorrange=ShStlims, strokecolor=(:black,1.0),strokewidth=1)
    end
    #Label(fig[2,1,Bottom()],"λ_"*string(n)*" = "*@sprintf("%.5E", evals[n]),fontsize = 32)

    #hidedecorations!(ax22)
    #hidespines!(ax22)

    colsize!(fig.layout,1,Aspect(1,1.0))

    colsize!(fig.layout,3,Aspect(1,1.0))


    Colorbar(fig[1,2],limits=colorrange=Plims,colormap=cgrad(:Blues_9, rev=true),flipaxis=true)
    Colorbar(fig[2,2],limits=colorrange=Tlims,colormap=:plasma,flipaxis=true)
    Colorbar(fig[1,4],limits=colorrange=Pefflims,colormap=:bwr,flipaxis=true)
    Colorbar(fig[2,4],limits=colorrange=ShStlims,colormap=:plasma,flipaxis=true)

    Label(fig[1,1,Bottom()],"Pressure",fontsize = 32, rotation=0)
    Label(fig[2,1,Bottom()],"Tension",fontsize = 32, rotation=0)

    Label(fig[1,3,Bottom()],"Effective Pressure",fontsize = 32, rotation=0)
    Label(fig[2,3,Bottom()],"Shear Stress",fontsize = 32, rotation=0)
    Label( fig[0,:],"Γ = "*string(params.γ)*", L₀ = "*string(params.L₀),fontsize = 32, color = (:black, 1))
    resize_to_layout!(fig)
    save(datadir(f,"eigenmodes","Stress.png"),fig)


end