using Catalyst
using DifferentialEquations
using LinearAlgebra
using CairoMakie

##
k₋ = 1/31e-3;
γ = .16;
tₘ = 50e-3;
λ(t) = γ*sin(2π*t/tₘ);
β = 1/26e-3;
k₊(t) = β*(1+λ(t));

rs = @reaction_network begin
    @variables k₊(t)
    (k₊,k₋), E <--> ES
end

params = (:k₋ => k₋)
u₀ = [:E => 1, :ES => 0]
tspan = (0,.8)

prob = DiscreteProblem(rs, u₀, tspan, params)
jump_prob = JumpProblem(rs, prob, Direct())
sol = solve(jump_prob, SSAStepper())

