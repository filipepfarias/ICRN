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
V   = 5e-18;                 # Original 1e-15
nₐ  = 6.023e24;              # Avogadro's number
k₁  = 1e6 / nₐ / V;          # 2nd order reaction
k₋₁ = 1e-4;                  # 1st order reaction 
k₂  = 0.1;                   # 1st order reaction 

K = [k₁;  # K₁
     k₋₁; # K₋₁
     k₂]; # K₂

# Initial conditions
ℰ, ℰ𝒜, 𝒜, ℬ = V * nₐ .* (5e-7, 0, 2e-7, 0);
ℰ, ℰ𝒜, 𝒜, ℬ = floor.(Int,(ℰ, ℰ𝒜, 𝒜, ℬ)) .+ 1 # Convertion to state-space index

𝛎 = Pr - Re;                  # Stoichiometric balance

𝗻ₖ = (30,30,30,30);           # State-space size

p₀ = ones(𝗻ₖ);
p₀ ./= sum(p₀); 
p₀[end] = 1 - sum(p₀[1:end-1]);
# p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;       # Initial condition

A = CMEOperator(𝛎,Re,K,𝗻ₖ);    # CME Operator