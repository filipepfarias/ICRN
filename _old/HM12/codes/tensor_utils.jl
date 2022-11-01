# ==== PACKAGES ====
using LinearAlgebra, MAT

# This line below just overloads a postfix operator to allow the element-wise transpose of a list of matrices
struct Transposer1 end; const ᴴ⁻ = Transposer1(); Base.:(*)(x, ::Transposer1) = [conj(y)' for y in x]
struct Transposer2 end; const ᵀ⁻ = Transposer2(); Base.:(*)(x, ::Transposer2) = [y' for y in x]
struct Transposer3 end; const ᴴ  = Transposer3(); Base.:(*)(x, ::Transposer3) = conj(x)'
struct Transposer4 end; const ᵀ  = Transposer4(); Base.:(*)(x, ::Transposer4) = x'

# ==== FUNCTIONS ====
Π(X...) = prod(X...);   # Alias for the product of a sequence
Σ(A...) = sum(A...)
function Base.:∘(a,b) reshape(hcat(([a].*b)...), (size(a)...,size(b)...)) end   # Outer product between two tensors

# The n-mode unfolding function
function unfold(𝓧, n)
    dims = size(𝓧); O = length(dims);   # Auxiliary variables
    nₓ  = (1:O)[1:O .≠ n]               # All dimensions but 'n'

    𝓧ₚ = permutedims(𝓧, [n, nₓ...])     # Re-order the fibers of 𝓧

    I₁ = dims[n]; I₂ = Π(dims[nₓ])      # Dimensions of unfolded matrix
    𝓧₍ₙ₎ = reshape(𝓧ₚ, (I₁, I₂))        # Unfolds the tensor into the matrix
end

# The n-mode folding function
function fold(𝓧₍ₙ₎, dims, n)
    O = length(dims);                   # The order of the tensor
    nₓ = [n, (1:O)[1:O.≠n]...]         # The unfolded dimensions ordering

    I = sortperm(nₓ)                    # Original dimensions ordering

    𝓧ₚ = reshape(𝓧₍ₙ₎, dims[nₓ])        # Reshape the matrix into a tensor
    𝓧  = permutedims(𝓧ₚ, I)             # Re-order the fibers of 𝓧
end

# The n-mode tensor product function 
function ×ₙ(𝓧, U; n=1)
    I = collect(size(𝓧)); I[n] = size(U)[1];  # Tensor product final dimensions
    I = tuple(I...)                           # Cast the array into a tuple                   

    𝓨ₙ = U * unfold(𝓧, n);                    # Performs the unfolded n-mode product
    𝓨  = fold(𝓨ₙ, I, n)                       # Folds the tensor to the original order
end

# The n-mode tensor product of a sequence of tensors
function ⨉ₙ(𝓧...; m=nothing)
    if m === nothing; m = 1:length(𝓧)-1; end; # Defines the order of n-mode products
    𝓨 = 𝓧[1];                                 # Initialize the 𝓨 tensor     
    for j in m                                # -
        𝓨 = ×ₙ(𝓨, 𝓧[1+j], n=j)                # Apply the tensor product in a sequence
    end                                       # -
    return 𝓨                                  # Return the resulting tensor
end

# The High-Order Singular Value Decomposition (HOSVD) + TruncatedHOSVD
function HOSVD(𝓧; R=nothing)
    if R === nothing; R = size(𝓧); end      # Defines the Rank-(R₁,⋯,Rₙ) to truncate

    U = Array[];                            # Creates an empty list of matrices
    for n = 1:ndims(𝓧)                      # - Iterate over each dimension -
        𝓧₍ₙ₎ = unfold(𝓧, n);                # Obtains unfolded matrix [𝓧]₍ₙ₎
        (U⁽ⁿ⁾,Σ⁽ⁿ⁾,V⁽ⁿ⁾) = svd(𝓧₍ₙ₎);       # Computes the SVD [𝓧]₍ₙ₎ = U⁽ⁿ⁾Σ⁽ⁿ⁾V⁽ⁿ⁾
        push!(U, U⁽ⁿ⁾);           # Adds first Rₙ columns of U⁽ⁿ⁾ to the list
    end                                     # - -

    𝓢 = ⨉ₙ(𝓧, (U)ᴴ⁻...)                     # Computes 𝓢 = 𝓧 ×₁ U⁽¹⁾ᵀ ×₂ ⋯ ×ₙ U⁽ᴺ⁾ᵀ
    return (𝓢, U)                           # Returns the decomposition
end

# The Higher-Order Orthogonal Iteration + Truncated HOOI
function HOOI(𝓧; R=nothing, It=500)
    if R === nothing; R = size(𝓧); end          # Defines the Rank-(R₁,⋯,Rₙ) to truncate
    N = ndims(𝓧)                                # Retrives the order of the tensor

    (_, U) = HOSVD(𝓧; R=R);                     # Initiates (U⁽¹⁾, ⋯, U⁽ᴺ⁾) using HOSVD
    for k = 1:It, n = 1:N                       # - Iterate over updates and each dimension -
        𝓤ₙ = ⨉ₙ(𝓧, (U)ᴴ⁻...; m=[1:n-1; n+1:N])  # Updates 𝓤ₙ = 𝓧 ⋯ ×ₙ₋₁ U⁽ⁿ⁻¹⁾ᵀ ×ₙ₊₁ U⁽ⁿ⁺¹⁾ᵀ ⋯ ×ₙ U⁽ᴺ⁾ᵀ

        𝓤₍ₙ₎ = unfold(𝓤ₙ, n);                   # Obtains unfolded matrix [𝓧]₍ₙ₎
        (U⁽ⁿ⁾,Σ⁽ⁿ⁾,V⁽ⁿ⁾) = svd(𝓤₍ₙ₎);           # Computes the SVD [𝓧]₍ₙ₎ = U⁽ⁿ⁾Σ⁽ⁿ⁾V⁽ⁿ⁾
        U[n] = U⁽ⁿ⁾[:,1:R[n]];                  # Updates (first Rₙ columns) of U⁽ⁿ⁾
    end                                         # - -

    𝓢 = ⨉ₙ(𝓧, (U)ᴴ⁻...)                          # Computes 𝓢 = 𝓧 ×₁ U⁽¹⁾ ×₂ ⋯ ×ₙ U⁽ᴺ⁾
    return (𝓢, U)                               # Returns the decomposition
end

# The plain-vanilla Alternating Least Squares (ALS) algorithm
function ALS(𝓧, R; ϵ=1e-16, It=5000)
    (I,J,K) = size(𝓧);                                          # Auxiliary variable
    (𝓧₍₁₎, 𝓧₍₂₎, 𝓧₍₃₎) = (unfold(𝓧, j) for j = 1:3 )           # Computes all unfoldings

    A⁽²⁾ᵢ₋₁ = randn(J,R); A⁽³⁾ᵢ₋₁ = randn(K,R); e₍ᵢ₋₁₎ = Inf;   # Initial A⁽²⁾, A⁽³⁾ and e₍ᵢ₎
    for i = 1:It                                                # - Iterate until converge -
        # Update the factor terms                               #
        global A⁽¹⁾ᵢ = 𝓧₍₁₎ * pinv((A⁽³⁾ᵢ₋₁ ⋄ A⁽²⁾ᵢ₋₁)ᵀ)        # Updates the estimate A⁽¹⁾
        global A⁽²⁾ᵢ = 𝓧₍₂₎ * pinv((A⁽³⁾ᵢ₋₁ ⋄ A⁽¹⁾ᵢ  )ᵀ)        # Updates the estimate A⁽²⁾
        global A⁽³⁾ᵢ = 𝓧₍₃₎ * pinv((A⁽²⁾ᵢ   ⋄ A⁽¹⁾ᵢ  )ᵀ)        # Updates the estimate A⁽³⁾
                                                                #
        # Checks for termination conditions                     #
        e₍ᵢ₎ = norm(𝓧₍₁₎ - A⁽¹⁾ᵢ*(A⁽³⁾ᵢ ⋄ A⁽²⁾ᵢ)ᵀ)              # Computes error measure e₍ᵢ₎
        if e₍ᵢ₋₁₎-e₍ᵢ₎ < ϵ; break;                              # Check for termination or update
        else              ; A⁽²⁾ᵢ₋₁ = A⁽²⁾ᵢ; A⁽³⁾ᵢ₋₁ = A⁽³⁾ᵢ; e₍ᵢ₋₁₎ = e₍ᵢ₎;
        end
    end                                                         # - -

    return (A⁽¹⁾ᵢ, A⁽²⁾ᵢ, A⁽³⁾ᵢ)                    # Returns the decomposition
end

# A function to compute the Rank-(R₁,⋯,Rₙ) and singular values |𝓢ₙ₌ₖ|_F
function rank_(𝓧; ϵ=1e-15)
    (𝓢,_) = HOSVD(𝓧)
    𝓢ₙ₌ₖ(n,k) = 𝓢[fill(:,n-1)..., k, fill(:,ndims(𝓢)-n)...]

    R = [ rank(unfold(𝓧, n), rtol=ϵ) for n = 1:ndims(𝓧) ]
    Σ = [ [norm(𝓢ₙ₌ₖ(n,k)) for k = 1:size(𝓢,n)] for n = 1:ndims(𝓢) ]
    
    return (R,Σ)
end

# Generalization the Kronecker product between two tensors
function ⊗(A,B) 
    C = zeros(size(A).*size(B)); S = size(A); 
    for i in CartesianIndices(A), j in CartesianIndices(B)
        k = [i[t]+S[t]*(j[t]-1)  for t = 1:ndims(A)]
        C[k...] = A[i] * B[j];
    end
    C
end 