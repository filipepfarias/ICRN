# Michaelis-Menten mechanism extracted from
# Qian and Ge on (9.1).

# ℰ <--> ℰ𝒜

specie = ["E", "EA"];

#     ℰ  ℰ𝒜  
Re = [ 1  0;   # k₁  ℰ  → ℰ𝒜
       0  1;   # k₋₁ ℰ𝒜  → ℰ
       0  1;   # k₂  ℰ  → ℰ𝒜
       1  0;   # k₋₂ ℰ𝒜  → ℰ
]
#     ℰ  ℰ𝒜 
Pr = [ 0  1 ;  # k₁  ℰ  → ℰ𝒜
       1  0 ;  # k₋₁ ℰ𝒜  → ℰ
       1  0 ;  # k₂  ℰ  → ℰ𝒜
       0  1 ;  # k₋₂ ℰ𝒜  → ℰ
]
# Deterministic to stochastic rate conversion
# k1 = 1 × 10⁶, k2 = 1 × 10⁻⁴, k3 = 0.1
# From Wilkinson, Stochastic Modelling for
# System Biology
# V   = 1e-15;                 # Original 1e-15
V   = 9e-17;                 # 
nₐ  = 6.022e23;              # Avogadro's number
k₁  = (5e-7) * 1e6;          # 1st order reaction
k₋₁ = 1e-4;                  # 1st order reaction 
k₂  = 0.1;                   # 1st order reaction 
k₋₂  = 1e-8 * 0;                   # 1st order reaction 

K = [k₁;  # K₁
     k₋₁; # K₋₁
     k₂;
     k₋₂;
     ]; # K₂

# Initial conditions
ℰ, ℰ𝒜 = (1, 1);
# ℰ, ℰ𝒜 = V * nₐ .* (2e-7, 0);
# ℰ, ℰ𝒜 = floor.(Int,(ℰ, ℰ𝒜)) .+ 1; # Convertion to state-space index
# S₀ = [ℰ, ℰ𝒜]' .- 1;
# ℰ, ℰ𝒜 = (ℰ-5:ℰ+5, ℰ𝒜-0:ℰ𝒜+0)

𝛎 = Pr - Re;                  # Stoichiometric balance

n = maximum(maximum.((ℰ, ℰ𝒜)));
𝗻ₖ = (n,n);                # State-space size

T = 0.0:.5:100.0;

###
I = Vector{Any}();
A = Matrix(operatorα(𝛎,Re,K,(2,2),[[0,1],[0,1]])[2:3,2:3]);
for f in 0.0:.05:1.0
       p₀ = [f, 1-f];
       pf = copy(p₀);
       Imut = Vector{Any}();
       for iT in eachindex(T)[2:end]
              P = exp(A*(T[iT]));
              PP = diagm(p₀)*P;
              imut = sum(PP .* log.( PP ./((P*p₀) .* p₀')  ));
              push!(Imut,imut);
       end
       push!(I,Imut);
end

eₚ = [entropy_production(sparse(p),sparse(A)) for p in pf_log];
d𝕊 = [-sum(sparselog(sparse(p)) .* (A*p)) for p in pf_log];
Eᵢₙ = [energy_input_rate(sparse(pf_log[end]),sparse(p),A) for p in pf_log];
hₑₓ = [entropy_flow(sparse(p),sparse(A)) for p in pf_log];
F = [free_energy(sparse(p),sparse(pf_log[end])) for p in pf_log];