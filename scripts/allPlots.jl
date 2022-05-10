# Import Julia packages
using DrWatson
@quickactivate

includet("$(projectdir())/scripts/analysisFunctions/allCurlsDivs.jl")
includet("$(projectdir())/scripts/analysisFunctions/allEigenModesLf.jl")
includet("$(projectdir())/scripts/analysisFunctions/allEigenModesLt.jl")
includet("$(projectdir())/scripts/analysisFunctions/allEigenModesLv.jl")
includet("$(projectdir())/scripts/analysisFunctions/internalVertexCurls.jl")
includet("$(projectdir())/scripts/analysisFunctions/intersections.jl")
includet("$(projectdir())/scripts/analysisFunctions/eigenvalues.jl")
includet("$(projectdir())/scripts/analysisFunctions/forceNeighbourhood.jl")
includet("$(projectdir())/scripts/analysisFunctions/forceNetwork.jl")
includet("$(projectdir())/scripts/analysisFunctions/fullSystem.jl")
includet("$(projectdir())/scripts/analysisFunctions/laplacianTableauLf.jl")
includet("$(projectdir())/scripts/analysisFunctions/laplacianTableauLt.jl")
includet("$(projectdir())/scripts/analysisFunctions/laplacianTableauLv.jl")
includet("$(projectdir())/scripts/analysisFunctions/orthogonality.jl")
includet("$(projectdir())/scripts/analysisFunctions/capitalPsivPotential.jl")
includet("$(projectdir())/scripts/analysisFunctions/capitalPsivSpectrum.jl")
includet("$(projectdir())/scripts/analysisFunctions/capitalPsivPotentialLtDerivative.jl")
includet("$(projectdir())/scripts/analysisFunctions/psicPotential.jl")
includet("$(projectdir())/scripts/analysisFunctions/psicPotentialLfDerivative.jl")
includet("$(projectdir())/scripts/analysisFunctions/psicSpectrum.jl")
includet("$(projectdir())/scripts/analysisFunctions/psivPotential.jl")
includet("$(projectdir())/scripts/analysisFunctions/psivPotentialLtDerivative.jl")
includet("$(projectdir())/scripts/analysisFunctions/psivSpectrum.jl")
includet("$(projectdir())/scripts/analysisFunctions/vertexCouples.jl")

centralCell=1
show=0

#dataDirs = ["data/AlexPaperParameters/$x" for x in readdir("data/AlexPaperParameters/") if isdir("data/AlexPaperParameters/$x")]
dataDirs = ["/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/AlexPaperParameters/$x" for x in readdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/AlexPaperParameters/$x") if isdir("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/VertexModelFigures/AlexPaperParameters/$x")]
dataDirs = ["data/figure7/$x" for x in readdir("data/figure7/$x") if isdir("data/figure7/$x")]

for dataDirectory in dataDirs
    allCurlsDivs(dataDirectory, show)
    allEigenModesLf(dataDirectory)
    # allEigenModesLt(dataDirectory)
    allEigenModesLv(dataDirectory)
    # internalVertexCurls(dataDirectory,show)
    capitalPsivPotential(dataDirectory,show)
    capitalPsivSpectrum(dataDirectory,show)
    # capitalPsivPotentialLtDerivative(dataDirectory,show)
    # eigenvalues(dataDirectory, show)
    # forceNeighbourhood(dataDirectory, centralCell, show)
    # forceNetwork(dataDirectory, centralCell, show)
    # fullSystem(dataDirectory, centralCell, 1, 1, 1, 0, 1, 0, 1, 0, show)
    laplacianTableauLf(dataDirectory, show)
    # laplacianTableauLt(dataDirectory, show)
    laplacianTableauLv(dataDirectory, show)
    # orthogonality(dataDirectory)
    psicPotential(dataDirectory, show)
    # psicPotentialLfDerivative(dataDirectory, show)
    psicSpectrum(dataDirectory, show)
    psivPotential(dataDirectory, show)
    # psivPotentialLtDerivative(dataDirectory, show)
    psivSpectrum(dataDirectory, show)
    # vertexCouples(dataDirectory, show)
    # intersectionDivsCurls(dataDirectory)
end
