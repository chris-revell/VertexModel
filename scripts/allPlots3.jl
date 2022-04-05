# Import Julia packages
using DrWatson
@quickactivate

includet("$(projectdir())/notebooks/allCurlsDivs.jl")
includet("$(projectdir())/notebooks/allEigenModesLf.jl")
includet("$(projectdir())/notebooks/allEigenModesLt.jl")
includet("$(projectdir())/notebooks/allEigenModesLv.jl")
# includet("$(projectdir())/notebooks/boundaryMidpointLinks.jl")
# includet("$(projectdir())/notebooks/intersections.jl")
includet("$(projectdir())/notebooks/eigenvalues.jl")
includet("$(projectdir())/notebooks/forceNeighbourhood.jl")
includet("$(projectdir())/notebooks/forceNetwork.jl")
includet("$(projectdir())/notebooks/fullSystem.jl")
includet("$(projectdir())/notebooks/laplacianTableauLf.jl")
includet("$(projectdir())/notebooks/laplacianTableauLv.jl")
includet("$(projectdir())/notebooks/orthogonality.jl")
includet("$(projectdir())/notebooks/psicPotential.jl")
includet("$(projectdir())/notebooks/psicPotentialDerivative.jl")
includet("$(projectdir())/notebooks/psicSpectrum.jl")
includet("$(projectdir())/notebooks/psivPotential.jl")
includet("$(projectdir())/notebooks/psivPotentialDerivative.jl")
includet("$(projectdir())/notebooks/psivSpectrum.jl")
includet("$(projectdir())/notebooks/vertexCouples.jl")

# dataDirs = String[]
# dataDirectory = "data/sims/2022-02-28-19-30-22"
# push!(dataDirs,"data/sims/2022-02-28-19-30-22")
# # dataDirectory = "data/sims/2022-02-28-19-30-22"
# push!(dataDirs,"data/sims/2022-02-28-19-30-22")
# # dataDirectory = "data/sims/2022-03-09-19-16-26"
# push!(dataDirs,"data/sims/2022-03-09-19-16-26")
# # dataDirectory = "data/sims/2022-03-15-15-47-00"
# push!(dataDirs,"data/sims/2022-03-15-15-47-00")
# # dataDirectory = "data/sims/2022-03-15-18-59-50"
# push!(dataDirs,"data/sims/2022-03-15-18-59-50")
# # dataDirectory = "data/sims/2022-03-16-12-26-01"
# push!(dataDirs,"data/sims/2022-03-16-12-26-01")
# # dataDirectory = "data/sims/2022-03-16-16-02-03"
# push!(dataDirs,"data/sims/2022-03-16-16-02-03")

centralCell=1
show=0

dataDirs = ["data/sims/$x" for x in readdir("data/sims/") if isdir("data/sims/$x")]

for dataDirectory in dataDirs
    allCurlsDivs(dataDirectory, show)
    allEigenModesLf(dataDirectory)
    allEigenModesLt(dataDirectory)
    allEigenModesLv(dataDirectory)
    # capitalPsivPotential(dataDirectory,show)
    eigenvalues(dataDirectory, show)
    forceNeighbourhood(dataDirectory, centralCell, show)
    forceNetwork(dataDirectory, centralCell, show)
    fullSystem(dataDirectory, centralCell, 1, 1, 1, 0, 1, 0, 1, 0, show)
    laplacianTableauLf(dataDirectory, show)
    laplacianTableauLv(dataDirectory, show)
    orthogonality(dataDirectory)
    psicPotential(dataDirectory, show)
    psicPotentialDerivative(dataDirectory, show)
    psicSpectrum(dataDirectory, show)
    psivPotential(dataDirectory, show)
    psivPotentialDerivative(dataDirectory, show)
    psivSpectrum(dataDirectory, show)
    vertexCouples(dataDirectory, show)
end
