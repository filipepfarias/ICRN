# ==== PACKAGES ====
using LinearAlgebra, MAT
include((@__DIR__)*"/HM12/codes/tensor_utils.jl");

# ==== FUNCTIONS ====
# The Tensor-Train product function
function ×¹(G...)
    # Auxiliary variables with dimension/rank information
    IR  = [size(G[1]); [size(G[i])[2:end] for i=2:length(G)]...]
    I,R = ([IR[i][1] for i=1:length(G)], [IR[i][2] for i=1:length(G)-1])
    
    # Tensor-train product loop
    𝓧 = zeros(I...)
    for i in CartesianIndices(𝓧), r in CartesianIndices((R[1],R[2],R[3]))
        𝓧[i] += G[1][i[1],r[1]] * G[2][r[1],i[2],r[2]] * G[3][r[2],i[3],r[3]] * G[4][r[3],i[4]]
    end
    return 𝓧    # Returns the resulting tensor
end

# The Tensor-Train Singular Value Decomposition (TTSVD) algorithm
function TTSVD(𝓧)
    I = size(𝓧)                                     # Auxiliary variable
                                                    #
    X¹ = reshape(𝓧, (I[1],Π(I[2:4])))               # Unfolds the tensor [𝓧]₁ ∈ ℝ(I₁,I₂I₃I₄)
    U,Σ,V = svd(X¹);    Σ = Diagonal(Σ);            # Computes the SVD [𝓧]₁ = UΣVᴴ
    R₁ = sum(diag(Σ) .> 1e-6)                       # Gets the number of nonzero singular-values
    U₁ = U[:,1:R₁];  V₁ = (Σ*(V)ᴴ)[1:R₁,:]          # Truncates pair (U,V) given R₁
    G₁ = U₁                                         # Stores the 1st factor G₁ ∈ ℝ(I₁,R₁)
                                                    #
    X² = reshape(V₁, (R₁*I[2], Π(I[3:4])))          # Unfolds the matrix X² = [V₁]₁ ∈ ℝ(R₁I₂,I₃I₄)
    U,Σ,V = svd(X²);    Σ = Diagonal(Σ);            # Computes the SVD X² = UΣVᴴ
    R₂ = sum(diag(Σ) .> 1e-6)                       # Gets the number of nonzero singular-values
    U₂ = U[:,1:R₂];  V₂ = (Σ*(V)ᴴ)[1:R₂,:]          # Truncates pair (U,V) given R₂
    𝓖₂ = reshape(U₂, (R₁, I[2], R₂))                # Reshapes the 2nd factor 𝓖₂ ∈ ℝ(R₁,I₂,R₂)
                                                    # 
    X³ = reshape(V₂, (R₂*I[3], I[4]))               # Unfolds the matrix X³ = [V₂]₁ ∈ ℝ(R₂I₃,I₄)
    U,Σ,V = svd(X³);    Σ = Diagonal(Σ);            # Computes the SVD X³ = UΣVᴴ
    R₃ = sum(diag(Σ) .> 1e-6)                       # Gets the number of nonzero singular-values
    U₃ = U[:,1:R₃];  V₃ = (Σ*(V)ᴴ)[1:R₃,:]          # Truncates pair (U,V) given R₃
    𝓖₃ = reshape(U₃, (R₂, I[3], R₃))                # Reshapes the 3rd factor 𝓖₃ ∈ ℝ(R₂,I₃,R₃)
    G₄ = V₃[1:R₃,:]                                 # Stores the 4th factor G₄ ∈ ℝ(R₃,I₄)

    return (G₁, 𝓖₂, 𝓖₃, G₄)                         # Returns the decomposition factors
end

# ==== SCRIPT ====
# -- Case 1: Generating 𝓧 from factors (G₁, 𝓖₂, 𝓖₃, G₄) --
# Randomly generates the TT factors and computes the tensor 𝓧 = G₁ ×¹ 𝓖₂ ×¹ 𝓖₃ ×¹ G₄
R = [4 3 4];     # Ranks used for the decomposition
G₁ = randn(5,R[1]); 𝓖₂ = randn(R[1],5,R[2]); 𝓖₃ = randn(R[2],5,R[3]); G₄ = randn(R[3],5);
𝓧 = ×¹(G₁, 𝓖₂, 𝓖₃, G₄);

G₁_, 𝓖₂_, 𝓖₃_, G₄_ = TTSVD(𝓧);      # Estimates the factors (G₁, 𝓖₂, 𝓖₃, G₄)

@assert 𝓧 ≈ ×¹(G₁_, 𝓖₂_, 𝓖₃_, G₄_)
@show norm(𝓧 - ×¹(G₁_, 𝓖₂_, 𝓖₃_, G₄_))^2

# -- Case 2: Directly generating 𝓧 randomly --
𝓧 = randn(5, 5, 5, 5);

G₁_, 𝓖₂_, 𝓖₃_, G₄_ = TTSVD(𝓧);      # Estimates the factors (G₁, 𝓖₂, 𝓖₃, G₄)

@assert 𝓧 ≈ ×¹(G₁_, 𝓖₂_, 𝓖₃_, G₄_)
@show norm(𝓧 - ×¹(G₁_, 𝓖₂_, 𝓖₃_, G₄_))^2

# ==== ====
