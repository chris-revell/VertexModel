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
using ColorSchemes


@from "$(projectdir())/src/VertexModelContainers.jl" using VertexModelContainers
@from "$(projectdir())/src/Laplacians.jl" using Laplacians
@from "$(projectdir())/src/AnalysisFunctions.jl" using AnalysisFunctions
@from "$(projectdir())/src/CellProperties.jl" using CellProperties

folder="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims/new_energy/100_cells/relaxed"
files=Glob.glob("new_energy/100_cells/relaxed/systemData*.jld2","C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims")
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
    g=vcat(cellPressures, -cellTensions)
    
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

    nv=LinRange(1, 2*nVerts, 2*nVerts)
    nc=LinRange(2*nVerts-2*nCells+1, 2*nVerts, 2*nCells)

    fig = Figure()
    ax=Axis(fig[1, 1], xlabel="n", ylabel="log₁₀λₙ", title="Γ = "*string(params.γ)*", L₀ = "*string(params.L₀))

    scatter!(ax,nv, log10.(sort(abs.(evalH))), color=:black,markersize=4, label=L"\lambda_n,\quad \mathcal{H}")
    scatter!(ax,nv, log10.(sort(abs.(evalLv))), color=:blue,markersize=4, label=L"\lambda_n,\quad \mathcal{L}_v^G")
    scatter!(ax,nc, log10.((abs.(reverse(sNF.^2)))), color=:red,markersize=4, label="σ²ₙ, N")
    
    #vlines!(ax,2*nVerts+1-((2*nCells)-1), color=:red)
    fig[1, 2] = Legend(fig, ax, framevisible = false)
    save(datadir(folder, "spectra","compare_spectra_log_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)
   
    fig = Figure()
    set_theme!(figure_padding=5, backgroundcolor=(:white,1.0), font="Helvetica", fontsize=19)
    ax=Axis(fig[1, 1], xlabel="Mode number, n", ylabel="λₙ", yscale=log10, title="Γ = "*string(params.γ)*", L₀ = "*string(params.L₀))
    hidedecorations!(ax, grid=true, ticks=false, label=false,ticklabels = false)
    #vspan!(197.5, 396.5, color = (:grey, 0.3))
    scatter!(ax,nv[4:end], ((abs.(evalH)))[4:end], color=:black,markersize=5, label=L"\lambda_n,\, \mathcal{H}")

    #scatter!(ax,nv, (sort(abs.(evalLv))), color=ColorSchemes.seaborn_colorblind6[1],markersize=4, label=L"\lambda_n,\, \mathcal{L}_v^G")


    scatter!(ax,nv[4:end], ((abs.(evmapLv)))[4:end], color=ColorSchemes.seaborn_colorblind6[2],markersize=5)
    scatter!(ax,nv[4:end], ((abs.(evmapgX)))[4:end], color=ColorSchemes.seaborn_colorblind6[3],markersize=5)


    elem_1 = [MarkerElement(color = :black, marker = :circle, markersize = 15)]

    elem_2 = [MarkerElement(color = ColorSchemes.seaborn_colorblind6[2], marker = :circle, markersize = 15)]

    elem_3 = [MarkerElement(color =ColorSchemes.seaborn_colorblind6[3], marker = :circle, markersize = 15)]

    # # elem_4 = [MarkerElement(color =ColorSchemes.seaborn_colorblind6[3], marker = :circle, markersize = 12)]
    Legend(fig[1, 2],
        [elem_1, elem_2, elem_3],
        [L"\lambda_n,\quad \mathcal{H}","Material","Geometric"],
        patchsize = (35, 35), rowgap = 10, framevisible = false)
    resize_to_layout!(fig)



    #fig[1, 2] = Legend(fig, ax, framevisible = false)
    save(datadir(folder,"spectra","Mode_contributions_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)
    
    mkpath(datadir(folder,"cell_modes"))
    mode_dir=mkpath(datadir(folder,"cell_modes", "Γ_"*string(params.γ)*"_L0_"*string(params.L₀)))

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
        save(datadir(mode_dir,"eigenmodes_Y$(@sprintf("%03d", 2*nCells-(n-1))).png"),fig)
    
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
    save(datadir(folder,"cell_plots","Area_Perimeter_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)

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
    save(datadir(folder,"cell_plots","circularity_shape_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)

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
    save(datadir(folder,"cell_plots","stress_pressure_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)

    modeSign=zeros(2*nCells)
    mode_mag_plus=zeros(2*nCells)
    mode_mag_minus=zeros(2*nCells)
    mode_magL_plus=zeros(2*nCells)
    mode_magL_minus=zeros(2*nCells)
    mode_magA_plus=zeros(2*nCells)
    mode_magA_minus=zeros(2*nCells)
    for n=1:2*nCells
        a=evecLc[1:nCells,n].*evecLc[nCells+1:end,n]
        bA=evecLc[1:nCells,n]
        bL=evecLc[nCells+1:end,n]
        modeSign[n]=length(a[a.>0])
        mode_mag_plus[n]=sum(a[a.>0])
        mode_mag_minus[n]=sum(a[a.<0])
        mode_magL_plus[n]=sum(bL[a.>0])
        mode_magL_minus[n]=sum(bL[a.<0])
        mode_magA_plus[n]=sum(bA[a.>0])
        mode_magA_minus[n]=sum(bA[a.<0])
    end
    n=LinRange(1, 2*nCells, 2*nCells)
    fig = Figure()
    ax=Axis(fig[1, 1], xlabel="n", ylabel="Count", title="Number of cells with the same sign A and L modes")
    
    #scatter!(ax,n, mode_magL_plus, color=:red, label="Σ L modes, same sign", markersize=7)
    scatter!(ax,n, modeSign, color=:black, label="Σ L modes, opposite", markersize=7)
    # scatter!(ax,n, mode_magA_plus.+ mode_magL_plus, color=:red, label="Σ same sign", markersize=7)
    # scatter!(ax,n, mode_magA_minus.+ mode_magL_minus, color=:blue, label="Σ opposite sign", markersize=7)
    
    #vlines!(ax,2*nVerts-(2*nCells) +0.5, color=:red)
    vlines!(ax,nCells+0.5, color=:red)
    #fig[1, 2] = Legend(fig, ax, framevisible = false)
    save(datadir(folder,"cell_modes","mode_sign_Γ_"*string(params.γ)*"_L0_"*string(params.L₀)*".png"),fig)
    


end