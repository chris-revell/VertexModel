# Import system data
  conditionsDict    = load("$dataDirectory/dataFinal.jld2")
  @unpack nVerts,nCells,nEdges,pressureExternal,γ,λ,viscousTimeScale,realTimetMax,tMax,dt,outputInterval,preferredPerimeter,preferredArea,outputTotal,realCycleTime,t1Threshold = conditionsDict["params"]
  matricesDict = load("$dataDirectory/matricesFinal.jld2")
  @unpack A,Aᵀ,B,Bᵀ,B̄,C,R,F,edgeTangents,edgeMidpoints,cellPositions,ϵ,cellAreas,boundaryVertices,edgeLengths = matricesDict["matrices"]

  T = makeCellLinks(conditionsDict["params"],matricesDict["matrices"])

  edgeTrapezia = makeEdgeTrapezia(conditionsDict["params"],matricesDict["matrices"])
  trapeziumAreas = abs.(area.(edgeTrapezia))

  linkTriangles = makeLinkTriangles(conditionsDict["params"],matricesDict["matrices"])
  linkTriangleAreas = abs.(area.(linkTriangles))

  cellPolygons = makeCellPolygons(conditionsDict["params"],matricesDict["matrices"])

  Lₜ = makeLt(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas,trapeziumAreas)

  eigenvectors = (eigen(Matrix(Lₜ))).vectors
  eigenvalues = (eigen(Matrix(Lₜ))).values


  vertexDivs = -1.0.*calculateVertexDivs(conditionsDict["params"],matricesDict["matrices"],T,linkTriangleAreas)

  onesVec = ones(nVerts)
  E = Diagonal(linkTriangleAreas)

  ḡ = ((onesVec'*E*vertexDivs)/(onesVec'*E*ones(nVerts))).*onesVec
  ğ = vertexDivs.-ḡ
  ψ̆ = zeros(nVerts)
  eigenmodeAmplitudes = Float64[]
  for k=2:nVerts
      numerator = -eigenvectors[:,k]'*E*ğ
      denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
      ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
      push!(eigenmodeAmplitudes,(numerator/denominator))
  end

  divLims = (-maximum(abs.(wideTildeVertexDivs)),maximum(abs.(wideTildeVertexDivs)))

  fig = Figure(resolution=(1000,1000),fontsize = 24)
  ax3 = Axis(fig[2,1][1,1],aspect=DataAspect(),fontsize=32)
  hidedecorations!(ax1)
  hidespines!(ax1)
  for k=1:nVerts
      poly!(ax3,linkTriangles[k],color=[wideTildeVertexDivs[k]],colorrange=divLims,colormap=:bwr,strokewidth=1,strokecolor=(:black,0.25)) #:bwr
  end
  for i=1:nCells
      poly!(ax3,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=1) #:bwr
  end
  Colorbar(fig[2,1][1,2],limits=divLims,colormap=:bwr,flipaxis=false,align=:left)
