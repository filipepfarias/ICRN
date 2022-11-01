# mat = Observable(Matrix{Float64}(undef,8,3))

# fig,ax,hm = heatmap(mat)

# for t in 0:2:400
#     pₜ = exp(Matrix(A)*t)*p₀[:]
#     # pₜ ./= sum(pₜ); 
#     mat[] = [sum(sum(reshape(pₜ,8,8,8),dims=3),dims=2)[:] sum(sum(reshape(pₜ,8,8,8),dims=3),dims=1)[:] sum(sum(reshape(pₜ,8,8,8),dims=1),dims=2)[:]]'
#     # hm = mat;
#     # fig
#     sleep(.0001)
# end


mat = Observable(Matrix{Float64}(undef,d,𝗻ₖ[1]))

fig,ax,hm = heatmap(mat)

for t in 0:1:400
    pₜ = reshape(exp(Matrix(A)*t)*p₀[:],𝗻ₖ)
    # pₜ ./= sum(pₜ); 
    mat[] = [sum(pₜ[circshift([n,repeat([:],d-1)...,],dd-1)...,]) for n in 1:𝗻ₖ[1], dd in 1:d]'
    # hm = mat;
    # fig
    sleep(.001)
end