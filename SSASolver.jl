using Pkg;
Pkg.activate(".");
using CME
using Random, Dates, FileIO, JLD2

model_nm = "MichaelisMenten"
model = "reactions/"*model_nm*".jl";
include(model);

path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")*"_SSA"
mkdir(path)

evol_𝕊 = zeros(length(T),1);

marg, marg_labels = CMEMarginals(𝗻ₖ,p₀,specie);
𝕊 = CMEEntropy(p₀);
𝔼 = CMEMean(𝗻ₖ,marg);
𝕍ar =CMEVariance(𝗻ₖ,𝔼,marg);

flname = path*"/"*model_nm*"_statistics_t"*string(0);
jldsave(flname, specie=specie,
marg_labels=marg_labels, 
marg=marg, E=𝔼, Var=𝕍ar, S=𝕊, t=0, T=T)
 
for iT in eachindex(T)
    sim = 0; max_sim = 2500;
    pS = zeros(𝗻ₖ...);

    for _ in 1:max_sim
        # global pS
        𝒮 = [rand(ℰ),ℰ𝒜,rand(𝒜),ℬ]' .-1;
        t,S = Gillespie(K, 𝛎, Re, 𝒮, T[iT]);
        pS[(S .+ 1)...] += 1;
    end

    pS ./= sum(pS);
    
    marg, marg_labels = CMEMarginals(𝗻ₖ,pS,specie);
    𝕊 = CMEEntropy(pS);
    i_marg = [1; 1 .+ cumsum(length(𝗻ₖ):-1:2)]
    𝔼 = CMEMean(𝗻ₖ,marg[i_marg]);
    𝕍ar =CMEVariance(𝗻ₖ,𝔼,marg[i_marg]);

    flname = path*"/"*model_nm*"_statistics_t"*string(iT);
    jldsave(flname, specie=specie,
    marg_labels=marg_labels, 
    marg=marg, E=𝔼, Var=𝕍ar, S=𝕊, t=T[iT], T=T)

    # evol_𝕊[iT] = 𝕊;
end


# println("Saving plots...")
include("misc_plotting_SSA.jl")
# using CairoMakie

# fig = Figure(resolution = (1600,1600));
# ax = Array{Any}(undef,length(specie),length(specie));
# idm = 1;
# for i in eachindex(specie), j in eachindex(specie)
#     global idm
#     if j > i
#         ax[i,j] = Axis(fig[i,j])
#         heatmap!(ax[i,j],marg[idm]')
#         idm += 1;
#     elseif j == i
#         ax[i,j] = Axis(fig[i,j],title=specie[i])
#         barplot!(ax[i,j],0:(𝗻ₖ[i]-1),marg[idm])
#         ylims!(ax[i,j],[0 1])
#         xlims!(ax[i,j],[0 𝗻ₖ[i]]);
#         idm += 1;
#     end
# end
# fig

# it = (T[1:end-1] .<= t') .& (T[2:end] .>= t');