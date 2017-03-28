using PyDSTool
using Base.Test

tic()
@time @testset "Linear Solve" begin include("linear.jl") end
@time @testset "Calcium Direct" begin include("calcium.jl") end
@time @testset "Calcium Parameterized" begin include("pf_calcium.jl") end
toc()
