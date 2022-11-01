
using GLMakie
using LinearAlgebra, SparseArrays
using DifferentialEquations

function J(νi,n) # νi per reaction
    return νi > 0 ? sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] : sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
end

function 𝗝(ν,n)
   return reduce(kron,J.(ν,n))
end

W(𝓘,Re,m) = reduce(kron,spdiagm.(eachcol(binomial.(𝓘,Re[m,:]')  )))

function CMEOperator(𝝼,Re,K,𝗻ₖ)
    𝓘 = hcat((:).(1,𝗻ₖ)...,);
    return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*W(𝓘,Re,m) for m in eachindex(𝝼[:,1])]))
end

d = 10;
Re = [Matrix(I(d-1)) zeros(Int64,d-1,1); zeros(Int64,d-1,1) reverse(Matrix(I(d-1)),dims=2)];
Pr = reverse(Re,dims=1);
K = 1e-2 * ones(2*(d-1));
𝗻ₖ = Tuple(repeat([10],d));

𝝼 = Pr-Re;
A = CMEOperator(𝝼,Re,K,𝗻ₖ)

p₀ = zeros(𝗻ₖ);
p₀[𝗻ₖ[1]] = 1;

f(u,p,t) = A*u;

p₀ = zeros(𝗻ₖ);
p₀[𝗻ₖ[1]] = 1;
u0 = p₀[:];

t=500.0;
prob = ODEProblem(f,u0, (0.0,t));
sol = solve(prob, Tsit5(), reltol=1e-15, abstol=1e-15, saveat=5)

mat = Observable(Matrix{Float64}(undef,d,𝗻ₖ[1]))

fig,ax,hm = heatmap(mat)

for i in eachindex(sol.u)
    pₜ = reshape(sol.u[i],𝗻ₖ)
    # pₜ ./= sum(pₜ); 
    mat[] = [sum(pₜ[circshift([n,repeat([:],d-1)...,],dd-1)...,]) for n in 1:𝗻ₖ[1], dd in 1:d]'
    # hm = mat;
    # fig
    sleep(.15)
end

# begin


# pₜ = reshape(sol.u[end],𝗻ₖ);
# heatmap([sum(pₜ[circshift([n,repeat([:],d-1)...,],dd-1)...,]) for n in 1:𝗻ₖ[1], dd in 1:d]')
# end