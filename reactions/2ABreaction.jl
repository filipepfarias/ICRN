# Example of reaction network extracted from 
# Erban and Chapman, on (1.55) and (1.56).

# 2𝒜     → ∅
# 𝒜 + ℬ → ∅
# ∅      → 𝒜
# ∅      → ℬ

specie = ["A", "B"];

#     𝒜  ℬ
Re = [2  0;   # K₁
      1  1;   # K₂
      0  0;   # K₃
      0  0];  # K₄

#     𝒜  ℬ
Pr = [0  0;   # K₁
      0  0;   # K₂
      1  0;   # K₃
      0  1];  # K₄

K = [1e-3;  # K₁
     1e-2;  # K₂
     1.2;   # K₃
     1];    # K₄
𝛎 = Pr - Re; # Stoichiometric balance

𝗻ₖ = (25,25); # State-space size
A = CMEOperator(𝛎,Re,K,𝗻ₖ);

p₀ = zeros(𝗻ₖ);
p₀[1,1] = 1.0;
p₀ ./= sum(p₀);