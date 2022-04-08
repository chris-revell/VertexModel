# Import Julia packages
using DrWatson
@quickactivate

includet("$(projectdir())/notebooks/allCurlsDivs.jl")
includet("$(projectdir())/notebooks/allEigenModesLf.jl")
includet("$(projectdir())/notebooks/allEigenModesLt.jl")
includet("$(projectdir())/notebooks/allEigenModesLv.jl")
includet("$(projectdir())/notebooks/internalVertexCurls.jl")
includet("$(projectdir())/notebooks/intersections.jl")
includet("$(projectdir())/notebooks/eigenvalues.jl")
includet("$(projectdir())/notebooks/forceNeighbourhood.jl")
includet("$(projectdir())/notebooks/forceNetwork.jl")
includet("$(projectdir())/notebooks/fullSystem.jl")
includet("$(projectdir())/notebooks/laplacianTableauLf.jl")
includet("$(projectdir())/notebooks/laplacianTableauLt.jl")
includet("$(projectdir())/notebooks/laplacianTableauLv.jl")
includet("$(projectdir())/notebooks/orthogonality.jl")
includet("$(projectdir())/notebooks/capitalPsivPotential.jl")
includet("$(projectdir())/notebooks/capitalPsivSpectrum.jl")
includet("$(projectdir())/notebooks/capitalPsivPotentialLtDerivative.jl")
includet("$(projectdir())/notebooks/psicPotential.jl")
includet("$(projectdir())/notebooks/psicPotentialLfDerivative.jl")
includet("$(projectdir())/notebooks/psicSpectrum.jl")
includet("$(projectdir())/notebooks/psivPotential.jl")
includet("$(projectdir())/notebooks/psivPotentialLtDerivative.jl")
includet("$(projectdir())/notebooks/psivSpectrum.jl")
includet("$(projectdir())/notebooks/vertexCouples.jl")

centralCell=1
show=0

dataDirs = ["data/sims/$x" for x in readdir("data/sims/") if isdir("data/sims/$x")]

for dataDirectory in dataDirs
    # allCurlsDivs(dataDirectory, show)
    # allEigenModesLf(dataDirectory)
    # allEigenModesLt(dataDirectory)
    # allEigenModesLv(dataDirectory)
    internalVertexCurls(dataDirectory,show)
    # capitalPsivPotential(dataDirectory,show)
    # capitalPsivSpectrum(dataDirectory,show)
    # capitalPsivPotentialLtDerivative(dataDirectory,show)
    # eigenvalues(dataDirectory, show)
    # forceNeighbourhood(dataDirectory, centralCell, show)
    # forceNetwork(dataDirectory, centralCell, show)
    # fullSystem(dataDirectory, centralCell, 1, 1, 1, 0, 1, 0, 1, 0, show)
    # laplacianTableauLf(dataDirectory, show)
    # laplacianTableauLt(dataDirectory, show)
    # laplacianTableauLv(dataDirectory, show)
    # orthogonality(dataDirectory)
    # psicPotential(dataDirectory, show)
    # psicPotentialLfDerivative(dataDirectory, show)
    # psicSpectrum(dataDirectory, show)
    # psivPotential(dataDirectory, show)
    # psivPotentialLtDerivative(dataDirectory, show)
    # psivSpectrum(dataDirectory, show)
    # vertexCouples(dataDirectory, show)
    # intersectionDivsCurls(dataDirectory)
end
