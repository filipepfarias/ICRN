using GLMakie

using LinearAlgebra, SparseArrays

function J(νi,n) # νi per reaction
    return νi > 0 ? sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] : sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
end

function 𝗝(ν,n)
   return reduce(kron,J.(ν,n))
end

# w(𝘅,r) = prod(binomial.(𝘅,r'),dims=2);
# W(𝘅,K,Re,m) = spdiagm(K[m]*w(𝘅,Re[m,:])[:]);
W(𝓘,Re,m) = reduce(kron,spdiagm.(eachcol(binomial.(𝓘,Re[m,:]')  )))
# spdiagm(K[m]*w(𝘅,Re[m,:])[:]);

function CMEOperator(𝝼,Re,K,𝗻ₖ)
    # d = length(𝗻ₖ);
    # 𝘅 = getindex.(CartesianIndices(𝗻ₖ)[:],(1:d)') .- 1;
    𝓘 = hcat((:).(1,𝗻ₖ)...,);
    return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*W(𝓘,Re,m) for m in eachindex(𝝼[:,1])]))
end

# 𝗻ₖ = (5,5,5,5,5);
# d = 2;
# Re = [I(d+1); zeros(Int64,d,d)];
# Pr = [zeros(Int64,d,d); I(d)]; 
# K = ones(2d,1);


# Re = [1 0 0; 0 1 0; 0 1 0; 0 0 1];
# Pr = [0 1 0; 0 0 1; 1 0 0; 0 1 0];
# K = 1e-3 * [1;2;2;1];
# 𝗻ₖ = (5,5,5);

# Re = [1 0; 0 1];
# Pr = [0 1; 1 0];
# K = 1 * [1;1];
# 𝗻ₖ = (20,20);

# Re = [1; 0];
# Pr = [0; 1];
# K = 1 * [1;1];
# 𝗻ₖ = Tuple(20);

# Re = [1 0 0; 0 1 0; 0 0 1; 0 1 0];
# Pr = [0 1 0; 0 0 1; 0 1 0; 1 0 0];
# K = 1e-2 * [1;1;1;1];
# 𝗻ₖ = (8,8,8);

d = 7;
Re = [Matrix(I(d-1)) zeros(Int64,d-1,1); zeros(Int64,d-1,1) reverse(Matrix(I(d-1)),dims=2)]
Pr = reverse(Re,dims=1);
K = 1e-2 * ones(2*(d-1));
𝗻ₖ = Tuple(repeat([4],d));

𝝼 = Pr-Re;
A = CMEOperator(𝝼,Re,K,𝗻ₖ)

pₜ₋₁ = p₀ = zeros(𝗻ₖ);
# pₜ₋₁ = p₀ = 1/25^2 * ones(𝗻ₖ);
p₀[𝗻ₖ[1]] = 1;
# p₀[1,20] = .5;

# p₀ /= sum(p₀);

# begin
# pₜ = exp(Matrix(A)*2)*p₀[:]
# heatmap([sum(sum(reshape(pₜ,8,8,8),dims=3),dims=2)[:] sum(sum(reshape(pₜ,8,8,8),dims=3),dims=1)[:] sum(sum(reshape(pₜ,8,8,8),dims=1),dims=2)[:]]')
# end
# heatmap(reshape(pₜ,8,8,8)[5,:,:] ./ sum(pₜ) )



# P[1] = 1;

# for iₜ in 1:Tₙ
#     # pₜ₋₁ = P[:,iₜ-1]
#     pₜ = pₜ₋₁ + A * pₜ₋₁ * .000001
#     # pₜ = pₜ ./ sum(pₜ)
#     # P[:,iₜ] = pₜ;
#     pₜ₋₁ = pₜ;
# end

# barplot(pₜ₋₁[:])