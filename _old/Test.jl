using CairoMakie

using LinearAlgebra, SparseArrays

function J(νi,n) # νi per reaction
    return νi > 0 ? sparse(I,n+νi,n+νi)[νi+1:end,1:end-νi] : sparse(I,n-νi,n-νi)[1:(end+νi),1-νi:end]
end

function 𝗝(ν,n)
   return reduce(kron,J.(ν,n))
end

H(𝘅,r) = prod(binomial.(𝘅,r'),dims=2);
# H(𝘅,r) = prod(binomial.(𝘅,r') .* factorial.(r'),dims=2);
W(𝘅,K,Re,m) = spdiagm(K[m]*H(𝘅,Re[m,:])[:]);

function CMEOperator(𝝼,Re,K,𝗻ₖ)
    d = length(𝗻ₖ);
    𝘅 = getindex.(CartesianIndices(𝗻ₖ)[:],(1:d)') .- 1;
    return (sum([(𝗝(-𝝼[m,:],𝗻ₖ[1]) - I)*W(𝘅,K,Re,m) for m in eachindex(𝝼[:,1])]))
end

# 𝘅 = getindex.(CartesianIndices((3,3))[:],(1:2)') .- 1

Re = [2 0; 1 1; 0 0; 0 0];
Pr = [0 0; 0 0; 1 0; 0 1];
K = [1e-3; 1e-2; 1.2; 1];
𝗻ₖ = (6,6);

# Re = [3; 2; 1; 0];
# Pr = [2; 3; 0; 1]; 
# K = [2.5e-4; 0.18; 37.5; 2200];
# 𝗻ₖ = Tuple(600);

# Re = [1 0 0; 0 1 1; 0 0 0; 0 1 1];
# Pr = [0 1 0; 1 0 1; 0 0 1; 0 1 0]; 
# K = [0.1; 1; 100; 100];
# 𝗻ₖ = (20,20,20);

# Re = [1; 0];
# Pr = [0; 1]; 
# K = [0.1; 1];
# 𝗻ₖ = Tuple(20);

𝝼 = Pr-Re;
A = CMEOperator(𝝼,Re,K,𝗻ₖ);

Tₙ = 100000; 
P = 1/prod(𝗻ₖ) * ones(prod(𝗻ₖ),Tₙ);
P[1] = 1;

# Tₙ = 1e7; 
# pₜ₋₁ = 1/prod(𝗻ₖ) * ones(prod(𝗻ₖ),1);

# Tₙ = 1e7; 
# pₜ₋₁ = 1/prod(𝗻ₖ) * ones(prod(𝗻ₖ),1);

# Tₙ = 1e7; 
# pₜ₋₁ = 1/prod(𝗻ₖ) * ones(prod(𝗻ₖ),1);

# p∞=exp(Matrix(A)*1)*p₀[:]
# P = exp(Matrix(A)*1e2)*p₀

for iₜ in 2:Tₙ
    pₜ₋₁ = P[:,iₜ-1]
    pₜ = pₜ₋₁ + A * pₜ₋₁ * .0002
    # pₜ = pₜ ./ sum(pₜ)
    P[:,iₜ] = pₜ;
end

# barplot(reshape(p∞ / sum(p∞),𝗻ₖ))
# barplot(reshape(p∞ / sum(p∞),𝗻ₖ))
# heatmap(reshape(p∞ / sum(p∞),(30,30)))
# f = Figure();
# ax = Axis(f[1, 1], title = "Step = .0001")

# lines!(ax,sum(P,dims=1)[:], title="Step = .1")
# f
