using DrWatson
@quickactivate
using VertexModel
using Glob
using DelimitedFiles
using UnPack
using FromFile
using JLD2

#generated from seven config

# using VertexModel; vertexModel(initialSystem="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims\\grow_layer\\L₀=1.0_realTimetMax=216000.0_t1Threshold=0.05_γ=0.01_23-11-16-14-08-14\\frameData\\systemData100.jld2", subFolder="vary_Gamma_smooth",realTimetMax=10*86400.0, γ=0.01,L₀=1)

# using VertexModel; vertexModel(initialSystem="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims\\grow_layer\\L₀=1.0_realTimetMax=216000.0_t1Threshold=0.05_γ=0.1_23-11-16-14-01-01\\frameData\\systemData100.jld2", subFolder="vary_Gamma_smooth",realTimetMax=10*86400.0, γ=0.1,L₀=1)

# using VertexModel; vertexModel(initialSystem="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims\\grow_layer\\L₀=1.0_realTimetMax=216000.0_t1Threshold=0.05_γ=1.0_23-11-16-14-03-13\\frameData\\systemData100.jld2", subFolder="vary_Gamma_smooth",realTimetMax=10*86400.0, γ=1,L₀=1)

# using VertexModel; vertexModel(initialSystem="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims\\grow_layer\\L₀=1.0_realTimetMax=216000.0_t1Threshold=0.05_γ=10.0_23-11-16-14-04-02\\frameData\\systemData100.jld2", subFolder="vary_Gamma_smooth",realTimetMax=10*86400.0, γ=10,L₀=1)

# using VertexModel; vertexModel(initialSystem="C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims\\grow_layer\\L₀=1.0_realTimetMax=216000.0_t1Threshold=0.05_γ=20.0_23-11-16-14-24-21\\frameData\\systemData100.jld2", subFolder="vary_Gamma_smooth",realTimetMax=10*86400.0, γ=20,L₀=1)

#generated from seven_eq config


f=Glob.glob("vary_Gamma_100_cells\\L₀=1.0_realTimetMax=346000.0_t1Threshold=0.05_γ=0.1_23-11-17-19-20-03/*/systemData100.jld2","C:\\Users\\v35431nc\\Documents\\VM_code\\VertexModel\\data\\sims" )[1]
T_len=[0.01, 0.075, 0.1, 0.2, 0.5]
for t in T_len
    @unpack R, matrices, params = load(f)
    @unpack γ, L₀ = params
    using VertexModel; vertexModel(initialSystem=f, subFolder="T1_len_relax",realTimetMax=25*86400.0, γ=γ,L₀=L₀,t1Threshold=t)

end
