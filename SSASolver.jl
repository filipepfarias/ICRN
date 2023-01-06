using CME

model_nm = "MichaelisMenten"
model = "reactions/"*model_nm*".jl";
include(model);

pS = p₀;
p = p₀;
sim = 0; max_sim = 5000;

while sim < max_sim
    𝒮 = [rand(ℰ),ℰ𝒜,rand(𝒜),ℬ]' .-1;
    t,S = Gillespie(K, 𝛎, Re, 𝒮, T[35]);
    pS[(S .+ 1)...] = 1;
    p += pS;
    pS[(S .+ 1)...] = 0;

    sim += 1;
end

p = p ./ sum(p);

# A = CMEOperator(𝛎,Re,K,𝗻ₖ); 
marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊, Si, Se = CMEStatistics(p[:],A,𝗻ₖ,specie);

# heatmap(marg[10]')
barplot(marg[10])

# it = (T[1:end-1] .<= t') .& (T[2:end] .>= t');