# DEPRECATED !!!

# using Pkg
# Pkg.activate(".")
## Plotting
using GLMakie, CairoMakie, FileIO, JLD2

path = "outputs/gv4Yl_20221117"
model_nm = "MichaelisMenten"

# For Michaelis-Menten reaction network
GLMakie.activate!()
fig = Figure(resolution = (1000,1000));

mat = [j > i ? 
        Observable(Matrix{Float64}(undef,𝗻ₖ[i],𝗻ₖ[j])) : 
        (j == i ? 
            Observable(Vector{Float64}(undef,𝗻ₖ[i])) : nothing)
        for i in eachindex(specie), j in eachindex(specie)];

ax = Array{Any}(undef,length(specie),length(specie));

for i in eachindex(specie), j in eachindex(specie)
    if j > i
        ax[i,j] = Axis(fig[i,j])
        heatmap!(ax[i,j],mat[i,j])
    elseif j == i
        ax[i,j] = Axis(fig[i,j],title=specie[i])
        barplot!(ax[i,j],0:(𝗻ₖ[i]-1),mat[i,j])
        ylims!(ax[i,j],[0 1])
        xlims!(ax[i,j],[0 𝗻ₖ[i]]);
    end
end
fig

try 
    mkpath(path*"/plots")
catch
    nothing
end

record(fig, path*"/plots/"*model_nm*"_anim.mp4", eachindex(T);
        framerate = 4) do iT
    iT -= 1;
    flname = path*"/"*model_nm*"_t"*string(iT)*"_marg_";
    for i in eachindex(specie), j in eachindex(specie)
        if j > i
            flsuffix = specie[i]*"_x_"*specie[j];
            mat[i,j][] = jldopen(flname*flsuffix)["p"]
        elseif i == j
            flsuffix = specie[i];
            ind = collect(eachindex(specie))
            mat[i,j][] = jldopen(flname*flsuffix)["p"];
        end
    end 
end

CairoMakie.activate!()
fig2 = Figure(resolution = (300,300));

flname = path*"/"*model_nm*"_mean";
𝔼 = jldopen(flname)["E"];
T = jldopen(flname)["T"];

fig2, ax, sp = series(T,𝔼, labels=specie);
axislegend(ax);
save(path*"/plots/"*model_nm*"_mean_evol.pdf", fig2, pt_per_unit = 2)

fig3 = Figure(resolution = (300,300));

flname = path*"/"*model_nm*"_entropy";
𝕊 = jldopen(flname)["S"];
d𝕊dt = jldopen(flname)["dSdt"];

fig3, ax, sp = series(T,[𝕊;d𝕊dt; 0 diff(𝕊[:])'], labels=["Entropy","Entropy balance","Test"]); 
axislegend(ax);
save(path*"/plots/"*model_nm*"_entrop_evol.pdf", fig3, pt_per_unit = 2)
