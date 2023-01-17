using CME
using Distributed

model_nm = "MichaelisMenten"
model = "reactions/"*model_nm*".jl";
include(model);

evol_𝕊 = zeros(length(T),1);
 
for iT in eachindex(T)
    global A, marg, evol_𝕊
    pS = zeros(size(p₀)...);
    sim = 0; max_sim = 500;

    for _ in 1:max_sim
        𝒮 = [rand(ℰ),ℰ𝒜,rand(𝒜),ℬ]' .-1;
        t,S = Gillespie(K, 𝛎, Re, 𝒮, T[iT]);
        pS[(S .+ 1)...] += 1;
    end

    pS ./= sum(pS);
    
    marg_labels, marg, = CMEMarginals(𝗻ₖ,p,specie);

    evol_𝕊[iT] = 𝕊;
end

p = p ./ sum(p);

# A = CMEOperator(𝛎,Re,K,𝗻ₖ); 
marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊, Si, Se = CMEStatistics(p[:],A,𝗻ₖ,specie);

# heatmap(marg[10]')
barplot(marg[10])

# it = (T[1:end-1] .<= t') .& (T[2:end] .>= t');