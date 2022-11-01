using CairoMakie

using LinearAlgebra, SparseArrays

function J(νi,n) # νi per reaction
    return νi > 0 ? sparse(I,n+νi,n+νi)[νi+1:end,1:end-νi] : sparse(I,n-νi,n-νi)[1:(end+νi),1-νi:end]
end

function 𝗝(ν,n)
   return reduce(kron,J.(ν,n))
end

H(𝘅,r) = prod(binomial.(𝘅,r'),dims=2);
W(𝘅,K,Re,m) = spdiagm(K[m]*H(𝘅,Re[m,:])[:]);

function CMEOperator(𝝼,Re,K,𝗻ₖ)
    d = length(𝗻ₖ);
    𝘅 = getindex.(CartesianIndices(𝗻ₖ)[:],(1:d)') .- 1;
    return (sum([(𝗝(-𝝼[m,:],𝗻ₖ[1]) - I)*W(𝘅,K,Re,m) for m in eachindex(𝝼[:,1])]))
end

Re = [1 0 0; 0 1 1; 0 0 0; 0 1 1];
Pr = [0 1 0; 1 0 1; 0 0 1; 0 1 0]; 
K = [0.1; 1; 100; 100];
𝗻ₖ = (20,20,20);

𝝼 = Pr-Re;
A = CMEOperator(𝝼,Re,K,𝗻ₖ);

Tₙ = 1e7; 
pₜ₋₁ = 1/prod(𝗻ₖ) * ones(prod(𝗻ₖ),1);
# pₜ₋₁ = zeros(prod(𝗻ₖ),1);
# pₜ₋₁[1] = 1;

# P[1] = 1;

for iₜ in 1:Tₙ
    # pₜ₋₁ = P[:,iₜ-1]
    pₜ = pₜ₋₁ + A * pₜ₋₁ * .000001
    # pₜ = pₜ ./ sum(pₜ)
    # P[:,iₜ] = pₜ;
    pₜ₋₁ = pₜ;
end

barplot(pₜ₋₁[:])
# barplot(P[:,end])
