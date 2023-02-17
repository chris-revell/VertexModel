using DrWatson
using DifferentialEquations
using CairoMakie
using StaticArrays

LJ(r) = 1/r^12 - 1/r^6

function model!(du,u,p,t)
    for (n,i) in u
        for (m,j) in u[i+1:end]
            F = LJ(norm(i-j)).*(i-j)/abs(i-j)
            du[i] -= F
            du[i+j] += F
        end
    end
    return du
end

u0 = @SVector [rand(2) for _ in 1:10]

prob = ODEProblem(model!,u0,(0.0,10.0))

sol = solve(prob)