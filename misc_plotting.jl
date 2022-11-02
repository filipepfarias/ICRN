## Plotting
using GLMakie


## For Michaelis-Menten reaction network
fig = Figure(resolution = (1000,1000));

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

for t in eachindex(p)
    𝓅  = reshape(p[t],𝗻ₖ...,);
    𝓅ₙ = sum(𝓅);
    for i in 1:4, j in 1:4
        if j > i
            d = deleteat!(collect(1:4), [i j])
            mat[i,j][] = reshape(sum(𝓅,dims=d) ./ 𝓅ₙ ,𝗻ₖ[i],𝗻ₖ[j])'
        elseif i == j
            ind = collect(1:4)
            mat[i,j][] = sum(𝓅,dims=deleteat!(ind,i))[:] ./𝓅ₙ ;
        end
    end
    # mat[] = reshape(p[t],𝗻ₖ);
    sleep(.5)
end

fig = Figure(resolution = (600,600));

𝔼 = zeros(length(𝗻ₖ),size(p)...,);

for t in eachindex(p)
    𝓅  = reshape(p[t],𝗻ₖ...,);
    𝓅ₙ = sum(𝓅);
    𝔼[:,t] = [
        sum(collect(0:(𝗻ₖ[i]-1)) .* sum(𝓅,dims=deleteat!(collect(1:length(𝗻ₖ)),i))[:] ./ 𝓅ₙ )
        for i in 1:length(𝗻ₖ)]   
end

fig, ax, sp = series(𝔼, labels=specie);
axislegend(ax);
fig

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
