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

using GLMakie

fig = Figure(resolution = (1600,1600));
ax = Array{Any}(undef,length(specie),length(specie));
idm = 1;
for i in eachindex(specie), j in eachindex(specie)
    global idm
    if j > i
        ax[i,j] = Axis(fig[i,j])
        heatmap!(ax[i,j],marg[idm]')
        idm += 1;
    elseif j == i
        ax[i,j] = Axis(fig[i,j],title=specie[i])
        barplot!(ax[i,j],0:(𝗻ₖ[i]-1),marg[idm])
        ylims!(ax[i,j],[0 1])
        xlims!(ax[i,j],[0 𝗻ₖ[i]]);
        idm += 1;
    end
end
fig

# it = (T[1:end-1] .<= t') .& (T[2:end] .>= t');