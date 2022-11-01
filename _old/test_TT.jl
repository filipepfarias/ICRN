
using GLMakie
using LinearAlgebra, SparseArrays, PROPACK
using DifferentialEquations

Π(X...) = prod(X...);

struct Transposer3 end; const ᴴ  = Transposer3(); Base.:(*)(x, ::Transposer3) = conj(x)'

function spunfold(a, dims::NTuple{2,Int})
    throw_dmrsa(dims, len) =
        throw(DimensionMismatch("new dimensions $(dims) must be consistent with array size $len"))

    i,j,v = findnz(a);
    oldI = [ LinearIndices(((:).(1,size(a))...,))[ij...] for ij in zip(i,j) ]
    lin = getindex.(CartesianIndices(dims)[oldI],(1:2)');
    return sparse(lin[:,1],lin[:,2],v)
end

# The n-mode unfolding function
function unfold(𝓧, n)
    dims = size(𝓧); O = length(dims);   # Auxiliary variables
    nₓ  = (1:O)[1:O .≠ n]               # All dimensions but 'n'

    𝓧ₚ = permutedims(𝓧, [n, nₓ...])     # Re-order the fibers of 𝓧

    I₁ = dims[n]; I₂ = Π(dims[nₓ])      # Dimensions of unfolded matrix
    𝓧₍ₙ₎ = reshape(𝓧ₚ, (I₁, I₂))        # Unfolds the tensor into the matrix

    return 𝓧₍ₙ₎
end

# The n-mode folding function
function fold(𝓧₍ₙ₎, dims, n)
    O = length(dims);                   # The order of the tensor
    nₓ = [n, (1:O)[1:O.≠n]...]         # The unfolded dimensions ordering

    I = sortperm(nₓ)                    # Original dimensions ordering

    𝓧ₚ = reshape(𝓧₍ₙ₎, dims[nₓ])        # Reshape the matrix into a tensor
    𝓧  = permutedims(𝓧ₚ, I)             # Re-order the fibers of 𝓧
end

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

d = 4;
Re = [Matrix(I(d-1)) zeros(Int64,d-1,1); zeros(Int64,d-1,1) reverse(Matrix(I(d-1)),dims=2)];
Pr = reverse(Re,dims=1);
K = 1e-2 * ones(2*(d-1));
𝗻ₖ = Tuple(repeat([4],d));

𝝼 = Pr-Re;
A = CMEOperator(𝝼,Re,K,𝗻ₖ);

function mTTSVD(X, dims)
    cores = Array{AbstractArray,1}(undef, length(dims)-1);

    Vᵢ = X;
    
    for i in 1:(length(dims)-1)
        Xⁱ = spunfold(Vᵢ, (dims[i],Π(dims[i+1:end])))  # Auxiliary variable
        U,Σ,V = tsvd(Xⁱ);    Σ = Diagonal(Σ);        
        V = sparse(V);                                     
        Rᵢ = sum(diag(Σ) .> 1e-6);
        Uᵢ = U[:,1:Rᵢ];  Vᵢ = (Σ*(V)ᴴ)[1:Rᵢ,:];
        Gᵢ = Uᵢ
        cores[i] = Gᵢ;
    end
    return cores
end

p₀ = zeros(𝗻ₖ);
p₀[𝗻ₖ[1]] = 1;

f(u,p,t) = A*u;

p₀ = zeros(𝗻ₖ);
p₀[𝗻ₖ[1]] = 1;
u0 = p₀[:];

# t=500.0;
# prob = ODEProblem(f,u0, (0.0,t));
# sol = solve(prob, Tsit5(), reltol=1e-15, abstol=1e-15, saveat=5)
