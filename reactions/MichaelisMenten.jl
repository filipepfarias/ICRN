# Michaelis-Menten mechanism extracted from
# Qian and Ge on (9.1).

# ℰ + 𝒜  → ℰ𝒜
# ℰ𝒜     → ℰ + 𝒜
# ℰ𝒜     → ℰ + ℬ

#     ℰ  ℰ𝒜 𝒜  ℬ 
Re = [1  0  1  0;   # k₁
      0  1  0  0;   # k₋₁
      0  1  0  0];  # k₂

#     ℰ  ℰ𝒜 𝒜  ℬ 
Pr = [0  1  0  0;   # k₁
      1  0  1  0;   # k₋₁
      1  0  0  1];  # k₂

# Deterministic to stochastic rate conversion
# k1 = 1 × 10⁶, k2 = 1 × 10⁻⁴, k3 = 0.1
# From Wilkinson, Stochastic Modelling for
# System Biology
V   = 3.5e-16;                 # Original 1e-15
nₐ  = 6.022e23;              # Avogadro's number
k₁  = 1e6 / nₐ / V;          # 2nd order reaction
k₋₁ = 1e-4;                  # 1st order reaction 
k₂  = 0.1;                   # 1st order reaction 

K = [k₁;  # K₁
     k₋₁; # K₋₁
     k₂]; # K₂

# Initial conditions
ℰ, ℰ𝒜, 𝒜, ℬ = V * nₐ .* (2e-7, 0, 5e-7, 0);
ℰ, ℰ𝒜, 𝒜, ℬ = floor.(Int,(ℰ, ℰ𝒜, 𝒜, ℬ)) .+ 1; # Convertion to state-space index
ℰ, ℰ𝒜, 𝒜, ℬ = (ℰ-15:ℰ+15, ℰ𝒜, 𝒜-15:𝒜+15, ℬ)
# ℰ, ℰ𝒜, 𝒜, ℬ = (7-2:7+2, 1, 16-2:16+2, 1);

𝛎 = Pr - Re;                  # Stoichiometric balance

𝗻ₖ = (128,128,128,128);           # State-space size

p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
# p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
p₀ ./= sum(p₀);

# p₀ = ones(𝗻ₖ);              # Uniform distribution
# p₀ ./= sum(p₀); 
# p₀[end] = 1 - sum(p₀[1:end-1]);

A = CMEOperator(𝛎,Re,K,𝗻ₖ);   # CME Operator      
