# Example of reaction network extracted from 
# Erban and Chapman, on (1.55) and (1.56).

# 2𝒜     → ∅
# 𝒜 + ℬ → ∅
# ∅      → 𝒜
# ∅      → ℬ
# 2𝒞     → ∅
# 𝒞 + 𝒟 → ∅
# ∅      → 𝒞
# ∅      → 𝒟

specie = ["A", "B", "C", "D"];

#     𝒜  ℬ  𝒞  𝒟
Re = [2  0   0  0;   # K₁
      1  1   0  0;   # K₂
      0  0   0  0;   # K₃
      0  0   0  0;  # K₄
      0  0   2  0;   # 
      0  0   1  1;   # 
      0  0   0  0;   # 
      0  0   0  0;];  # 

#     𝒜  ℬ  𝒞  𝒟
Pr = [0  0   0  0;   # K₁
      0  0   0  0;   # K₂
      1  0   0  0;   # K₃
      0  1   0  0;  # K₄
      0  0   0  0;   # 
      0  0   0  0;   # 
      0  0   1  0;   # 
      0  0   0  1;]  # 

K = [1e-3;  # K₁
     1e-2;  # K₂
     1.2;   # K₃
     1;    # K₄
     1e-3;  # K₁
     1e-2;  # K₂
     1.2;   # K₃
     1];
𝛎 = Pr - Re; # Stoichiometric balance

𝗻ₖ = (16,16,16,16); # State-space size
A = CMEOperator(𝛎,Re,K,𝗻ₖ);

p₀ = zeros(𝗻ₖ);
p₀[1,1,1,1] = 1.0;
p₀ ./= sum(p₀);