#
#  Parameters.jl
#  VertexModelJL
#
#  Created by Christopher Revell on 15/02/2021.
#
#
#

module Parameters

gamma               = 0.172                 # Parameters in energy relaxation. Hard wired from data.    (0.05#0.172 ??)
lamda               = -0.259                # Parameters in energy relaxation. Hard wired from data.    (-0.1#-0.259 ??)
tStar               = 20.0                  # Relaxation rate. Approx from Sarah's data.
realTimetMax        = 100.0                 # Real time maximum system run time /seconds
tMax                = realTimetMax/tStar    # Non dimensionalised maximum system run time
dt                  = 0.01                  # Non dimensionalised time step
outputInterval      = tMax/100.0            # Time interval for storing system data (non dimensionalised)
meanCellCycle       = 4.0*(60^2)            # Real time cell cycle time /seconds
nonDimCellCycle     = meanCellCycle/tStar   # Non dimensionalised cell cycle time? (Parameter used in cell division time)
preferredPerimeter  = -lamda/(2*gamma)      # Cell preferred perimeter

export gamma,lamda,tStar,realTimetMax,tMax,dt,outputInterval,meanCellCycle,nonDimCellCycle,alphaPrefA,preferredPerimeter

end
