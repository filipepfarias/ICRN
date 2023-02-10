# Michaelis-Menten mechanism extracted from
# Qian and Ge on (9.1).

# ℰ + 𝒜  → ℰ𝒜
# ℰ𝒜     → ℰ + 𝒜
# ℰ𝒜     → ℰ + ℬ

specie = ["E", "EA", "A", "B"];

#     ℰ  ℰ𝒜 𝒜  ℬ 
Re = [1  0  1  0;   # k₁  ℰ + 𝒜  → ℰ𝒜
      0  1  0  0;   # k₋₁ ℰ𝒜  → ℰ + 𝒜
      0  1  0  0;   # k₂  ℰ𝒜  → ℰ + ℬ
      1  0  0  1;   # k₋₂ ℰ + ℬ →  ℰ𝒜
      0  0  0  0;  # k₃
      1  0  0  0;  # k₋₃
      0  0  0  0;  # k₄
      0  0  1  0;  # k₋₄
      0  0  0  0;  # k₅
      0  1  0  0;  # k₋₅
      0  0  0  0;  # k₆
      0  0  0  1;  # k₋₆
]
#     ℰ  ℰ𝒜 𝒜  ℬ 
Pr = [0  1  0  0;   # k₁  ℰ + 𝒜  → ℰ𝒜
      1  0  1  0;   # k₋₁ ℰ𝒜  → ℰ + 𝒜
      1  0  0  1;  # k₂   ℰ𝒜  → ℰ + ℬ
      0  1  0  0;  # k₋₂  ℰ + ℬ →  ℰ𝒜
      1  0  0  0;  # k₃
      0  0  0  0;  # k₋₃
      0  0  1  0;  # k₄
      0  0  0  0;  # k₋₄
      0  1  0  0;  # k₅
      0  0  0  0;  # k₋₅
      0  0  0  1;  # k₆
      0  0  0  0;  # k₋₆
]
# Deterministic to stochastic rate conversion
# k1 = 1 × 10⁶, k2 = 1 × 10⁻⁴, k3 = 0.1
# From Wilkinson, Stochastic Modelling for
# System Biology
# V   = 1e-15;                 # Original 1e-15
V   = 5e-17;                 # 
nₐ  = 6.022e23;              # Avogadro's number
k₁  = 1e6 / nₐ / V;          # 2nd order reaction
k₋₁ = 1e-4;                  # 1st order reaction 
k₂  = 0.1;                   # 1st order reaction 
k₋₂  = 0.0001;                   # 1st order reaction 
k₃ = 0.0;
k₋₃ = 0.0;
k₄ = 0.0;
k₋₄ = 0.0;
k₅ = 0.0;
k₋₅ = 0.0;
k₆ = 0.0;
k₋₆ = 0.0;

K = [k₁;  # K₁
     k₋₁; # K₋₁
     k₂;
     k₋₂;
     k₃
     k₋₃;
     k₄;
     k₋₄;
     k₅;
     k₋₅;
     k₆;
     k₋₆;
     ]; # K₂

# Initial conditions
ℰ, ℰ𝒜, 𝒜, ℬ = V * nₐ .* (2e-7, 0, 5e-7, 0);
ℰ, ℰ𝒜, 𝒜, ℬ = floor.(Int,(ℰ, ℰ𝒜, 𝒜, ℬ)) .+ 1; # Convertion to state-space index
𝒮 = (ℰ-5:ℰ+5, ℰ𝒜-0:ℰ𝒜+0, 𝒜-5:𝒜+5, ℬ-0:ℬ+0);
S₀ = [ℰ, ℰ𝒜, 𝒜, ℬ]' .- 1;

𝛎 = Pr - Re;                  # Stoichiometric balance
n = maximum(maximum.(𝒮));
𝗻ₖ = (Int(floor(3*n/5)),Int(floor(3*n/5)),n,n);                # State-space size

T = 0.0:.5:100.0;