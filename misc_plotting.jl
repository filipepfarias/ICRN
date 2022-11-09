## Plotting
using GLMakie, CairoMakie, FileIO, JLD2

# path = "outputs/4x3Po_20221109"

## For Michaelis-Menten reaction network
GLMakie.activate!()
fig = Figure(resolution = (1600,1600));

specie = ["E", "EA", "A", "B"];

mat = [j > i ? 
        Observable(Matrix{Float64}(undef,𝗻ₖ[i],𝗻ₖ[j])) : 
        (j == i ? 
            Observable(Vector{Float64}(undef,𝗻ₖ[i])) : nothing)
        for i in 1:4, j in 1:4];

ax = Array{Any}(undef,4,4);

for i in 1:4, j in 1:4
    if j > i
        ax[i,j] = Axis(fig[i,j])
        heatmap!(ax[i,j],mat[i,j])
    elseif j == i
        ax[i,j] = Axis(fig[i,j],title=specie[i])
        barplot!(ax[i,j],0:(𝗻ₖ[i]-1),mat[i,j])
        ylims!(ax[i,j],[0 1])
    end
end
fig

try 
    mkdir(path*"/plots")
    catch
        nothing
    end

record(fig, path*"/plots/MichaelisMenten_anim.mp4", eachindex(T);
        framerate = 4) do iT
    iT -= 1;
    flname = path*"/MichaelisMenten_t"*string(iT)*"_marg_";
    for i in 1:4, j in 1:4
        if j > i
            flsuffix = specie[i]*"_x_"*specie[j];
            mat[i,j][] = jldopen(flname*flsuffix)["p"]
        elseif i == j
            flsuffix = specie[i];
            ind = collect(1:4)
            mat[i,j][] = jldopen(flname*flsuffix)["p"];
        end
    end 
end

CairoMakie.activate!()
fig2 = Figure(resolution = (300,300));

flname = path*"/MichaelisMenten_mean";
𝔼 = jldopen(flname)["E"];

fig2, ax, sp = series(𝔼, labels=specie);
axislegend(ax);
save(path*"/plots/MichaelisMenten_mean_evol.pdf", fig2, pt_per_unit = 2)

# # Entropy
# g = Figure();
# ax = Axis(g[1,1], yscale=log10)
# ylims!(ax, 0.1,100) 
# lines!(ax,)
# current_figure()

# g = Figure();
# ax = Axis(g[1,1], yscale=log10)
# ylims!(ax, 0.01,10) 
# lines!(ax,)
# current_figure()
