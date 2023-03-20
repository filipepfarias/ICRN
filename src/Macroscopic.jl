# using DifferentialEquations: solve, ODEProblem, RK4

function φₘ(𝐗, Re, kₘ, V = 1)
    # Qian, Ge (2021) Eq. 5.20
    return kₘ*V*prod(α(𝐗,Re)./(V.^Re))
end

# function φ̂ₘ(𝐱, Re, kₘ)
#     # Qian, Ge (2021) Eq. 5.20
#     return kₘ*prod((^).(𝐱,Re))
# end

function R(𝐱,m)  # Concentration Flux
    global K, Re
    k = rate_conversion(K[m],Re[m,:]);
    return k*prod((^).(𝐱,Re[m,:]))
end

# function d𝐱(𝐱,p,t)
#     global Re, 𝛎, K
#     J = [φ̂ₘ(𝐱,Re[:,m],K[m]) for m in 1:size(K,1)];
#     return sum(J[1:2:end]-J[2:2:end])[:]
# end

function d𝐱dt!(𝐱,p,t)
    global K
    return sum([𝛎[m,n]*R(𝐱,m) for n in 1:size(𝛎,2), m in 1:size(K,1)],dims=2)[:]
end

function macroscopic_entropy_production(𝐱)
    Rₘ = [(R(𝐱,m)-R(𝐱,2m)) for m in 1:Int(size(K,1)/2)]';
    logRₘ = [((R(𝐱,m) > 0)*log(R(𝐱,m)) - (R(𝐱,2m) > 0)*log(R(𝐱,2m))) for m in 1:Int(size(K,1)/2)];
    return sum(Rₘ * logRₘ);
end

function generalized_free_energy(𝐱,𝐱ₛₛ)
    Rₘ = -[(R(𝐱,m)-R(𝐱,2m)) for m in 1:Int(size(K,1)/2)]';
    logRₘ = -[((R(𝐱,m) > 0)*log(R(𝐱,m)) - (R(𝐱,2m) > 0)*log(R(𝐱,2m))) for m in 1:Int(size(K,1)/2)];
    logRₛₛ = [((R(𝐱ₛₛ,m) > 0)*log(R(𝐱ₛₛ,m)) - (R(𝐱ₛₛ,2m) > 0)*log(R(𝐱ₛₛ,2m))) for m in 1:Int(size(K,1)/2)];
    return sum(Rₘ * (logRₘ + logRₛₛ));
end

function macroscopic_energy_input(𝐱,𝐱ₛₛ)
    Rₘ = -[(R(𝐱,m)-R(𝐱,2m)) for m in 1:Int(size(K,1)/2)]';
    logRₘ = [((R(𝐱,m) > 0)*log(R(𝐱,m)) - (R(𝐱,2m) > 0)*log(R(𝐱,2m))) for m in 1:Int(size(K,1)/2)];
    logRₛₛ = -[((R(𝐱ₛₛ,m) > 0)*log(R(𝐱ₛₛ,m)) - (R(𝐱ₛₛ,2m) > 0)*log(R(𝐱ₛₛ,2m))) for m in 1:Int(size(K,1)/2)];
    return sum(Rₘ * (logRₘ + logRₛₛ));
end

function rate_conversion(k,Reₘ)
    order = sum(Reₘ);
    s = prod(Reₘ[Reₘ .!= 0])
    v = -min(1,order-1);
    return s\k*(nₐ*V)^-v
end

# function DetSolver(path, model_nm; molecules=false)
#     global specie, Κ
#     model = "reactions/"*model_nm*".jl";
#     include(model);

#     p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
#     p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
#     # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
#     p₀ ./= sum(p₀);

#     # p₀ = ones(𝗻ₖ);              # Uniform distribution
#     # p₀ ./= sum(p₀); 
#     # p₀[end] = 1 - sum(p₀[1:end-1]);

#     Κ = K;
#     Κ[1] = K[1] * nₐ * V; #Molecules to concentrations conversion (Wilkinson)

#     # function g!(𝘅,p,t)
#     #     d𝘅[1] = (Κ[2]+Κ[3]) * 𝘅[2] - Κ[1]*𝘅[1]*𝘅[3]
#     #     d𝘅[2] = Κ[1]*𝘅[1]*𝘅[3] - (Κ[2]+Κ[3]) * 𝘅[2]
#     #     d𝘅[3] = Κ[2] * 𝘅[2] - Κ[1]*𝘅[1]*𝘅[3]
#     #     d𝘅[4] = Κ[3] * 𝘅[2] 
#     # end

#     function f!(𝘅,p,t)
#         global specie, Κ
#         d𝘅 = sum([𝛎[ℓ,i]*Jflux(𝘅,ℓ) for i in eachindex(specie), ℓ in eachindex(Κ)],dims=2)[:]
#     end

#     𝘅₀ = [2e-7; 0; 5e-7; 0]

#     prob = ODEProblem(f!,𝘅₀,(T[1],T[end]))
#     sol = solve(prob,RK4();dt= .5,adaptive=false)

#     x = molecules ? hcat(sol.u...) * nₐ * V : hcat(sol.u...);

#     mkpath(path)
#     flname = path*"/"*model_nm;
#     jldsave(flname, specie=specie, x=x, T=T)

#     return (specie, x, T)
# end



function GibbsFreeEnergy(𝘅)    # Gibbs free energy
    g1 = [(-1)^ℓ * Jflux(𝘅,ℓ) for ℓ in eachindex(K)]
    g2 = [-(-1)^ℓ * log(Jflux(𝘅,ℓ)) for ℓ in eachindex(K)]
    g1[isinf.(g1) .|| g1 .== 0.0 ] .= 0.0;
    g2[isinf.(g2) .|| g2 .== 0.0 ] .= 0.0;
    g1g2 = sum(g1 * g2')
    d𝘅 = sum(g1g2[!isnan.(g1g2)]);
    return d𝘅
end