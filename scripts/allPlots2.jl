# Import Julia packages
using DrWatson
@quickactivate

# notebooks = [f for f in readdir("$(projectdir())/notebooks") if f[end-2:end]==".jl"]

# dataDirectory = "data/sims/2022-02-28-19-30-22"
# dataDirectory = "data/sims/2022-02-28-19-30-22"
# dataDirectory = "data/sims/2022-03-09-19-16-26"
# dataDirectory = "data/sims/2022-03-15-15-47-00"
# dataDirectory = "data/sims/2022-03-15-18-59-50"
# dataDirectory = "data/sims/2022-03-16-12-26-01"
dataDirectory = "data/sims/2022-03-16-16-02-03"


include("$(projectdir())/notebooks/allCurlsDivs.jl")
include("$(projectdir())/notebooks/allEigenModesLf.jl")
include("$(projectdir())/notebooks/allEigenModesLv.jl")
# include("$(projectdir())/notebooks/boundaryMidpointLinks.jl")
include("$(projectdir())/notebooks/eigenvalues.jl")
include("$(projectdir())/notebooks/forceNeighbourhoodOnly.jl")
include("$(projectdir())/notebooks/forceNetworkOnly.jl")
include("$(projectdir())/notebooks/fullSystem.jl")
include("$(projectdir())/notebooks/laplacianTableauLf.jl")
include("$(projectdir())/notebooks/laplacianTableauLv.jl")
# include("$(projectdir())/notebooks/orthogonality.jl")
include("$(projectdir())/notebooks/psicPotential.jl")
include("$(projectdir())/notebooks/psicPotentialLfDerivative.jl")
include("$(projectdir())/notebooks/psicSpectrum.jl")
include("$(projectdir())/notebooks/psivPotential.jl")
include("$(projectdir())/notebooks/psivPotentialLtDerivative.jl")
include("$(projectdir())/notebooks/psivSpectrum.jl")
include("$(projectdir())/notebooks/vertexCouples.jl")
