# Example of reaction network extracted from 
# Erban and Chapman, on (1.55) and (1.56).

# 3𝒜 → 2𝒜
# 2𝒜 → 3𝒜
# 𝒜  → ∅
# ∅  → 𝒜

specie = ["A"];

#     𝒜
Re = [3;   # K₁
      2;   # K₂
      1;   # K₃
      0];  # K₄

#     𝒜
Pr = [2;   # K₁
      3;   # K₂
      0;   # K₃
      1];  # K₄

K = [2.5e-4;   # K₁
     .18;      # K₂
     37.5;     # K₃
     2200];    # K₄
𝛎 = Pr - Re; # Stoichiometric balance

𝗻ₖ = Tuple(600); # State-space size
A = CMEOperator(𝛎,Re,K,𝗻ₖ);

p₀ = ones(𝗻ₖ);
# p₀[1] = 1.0;
p₀ ./= sum(p₀);