### VertexModel.jl

Open Julia REPL at project root directory.

Activate environment with `using Pkg; Pkg.activate(".")` or access package manager with `]` and then enter `activate .` before exiting package manager with backspace. Alternatively, load DrWatson library with `using DrWatson` and enter `@quickactivate`.

If running for the first time, run `Using Pkg; Pkg.instantiate()` or access package manager with `]` and then enter `instantiate` before exiting package manager with backspace.

Load test parameters by entering `includet("scripts/testParameters.jl"). This instantiates variables with values for running the vertex model code. 

These steps can be made easier in future by running `include("loadAll.jl")`.

Once `VertexModel.jl` has been precompiled, run function `vertexModel` using test parameters by entering:
`vertexModel(initialSystem,realTimetMax,realCycleTime,γ,L₀,A₀,viscousTimeScale,dt,pressureExternal,t1Threshold,outputTotal,outputToggle,plotToggle;subFolder="test")`

TODO:

- Control distribution of cell divisions by selecting from random distribution.
- How do we choose T1 transition edge length threshold? Should this depend on preferred cell perimeter? 
- Remove vertices where there are two short peripheral edges belonging to the same cell.
- Consolidate ordering of vertices and edges around cells using topology into a function that can be used everywhere; reduce reliance on calculating angles. 
