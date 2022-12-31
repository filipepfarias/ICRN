using CME

model_nm = "MichaelisMenten"
model = "reactions/"*model_nm*".jl";
include(model);

p₀ = [];

t,S = Gillespie(K, 𝛎, S₀, T[end])