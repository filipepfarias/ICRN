using CME
using DifferentialEquations

model_nm = "MichaelisMenten"
model = "reactions/"*model_nm*".jl";
include(model);
Κ = K;
Κ[1] = K[1] * nₐ * V; #Molecules to concentrations conversion (Wilkinson)
𝒥(𝘅,ℓ) = Κ[ℓ]*prod((^).(𝘅,Re[ℓ,:]))

# function g!(𝘅,p,t)
#     d𝘅[1] = (Κ[2]+Κ[3]) * 𝘅[2] - Κ[1]*𝘅[1]*𝘅[3]
#     d𝘅[2] = Κ[1]*𝘅[1]*𝘅[3] - (Κ[2]+Κ[3]) * 𝘅[2]
#     d𝘅[3] = Κ[2] * 𝘅[2] - Κ[1]*𝘅[1]*𝘅[3]
#     d𝘅[4] = Κ[3] * 𝘅[2] 
# end

function f!(𝘅,p,t)
    d𝘅 = sum([𝛎[ℓ,i]*𝒥(𝘅,ℓ) for i in eachindex(specie), ℓ in eachindex(K)],dims=2)[:]
end

𝘅₀ = [2e-7; 0; 5e-7; 0]

tspan = (0.0,100.0)
prob = ODEProblem(f!,𝘅₀,tspan)
sol = solve(prob,RK4();dt= .5,adaptive=false)

function G(𝘅)
    g1 = [(-1)^ℓ * 𝒥(𝘅,ℓ) for ℓ in eachindex(K)]
    g2 = [-(-1)^ℓ * log(𝒥(𝘅,ℓ)) for ℓ in eachindex(K)]
    d𝘅 = sum(g1 * g2')
    return d𝘅
end

tspan = (0.0,100.0)
prob = ODEProblem(G!,𝘅₀,tspan)
sol = solve(prob,RK4();dt= .5/10,adaptive=false)

using GLMakie