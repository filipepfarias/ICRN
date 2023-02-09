using DifferentialEquations: solve, ODEProblem, RK4

function DetSolver(path, model_nm; molecules=false)
    global specie, Κ
    model = "reactions/"*model_nm*".jl";
    include(model);

    p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
    p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
    # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
    p₀ ./= sum(p₀);

    # p₀ = ones(𝗻ₖ);              # Uniform distribution
    # p₀ ./= sum(p₀); 
    # p₀[end] = 1 - sum(p₀[1:end-1]);

    Κ = K;
    Κ[1] = K[1] * nₐ * V; #Molecules to concentrations conversion (Wilkinson)

    # function g!(𝘅,p,t)
    #     d𝘅[1] = (Κ[2]+Κ[3]) * 𝘅[2] - Κ[1]*𝘅[1]*𝘅[3]
    #     d𝘅[2] = Κ[1]*𝘅[1]*𝘅[3] - (Κ[2]+Κ[3]) * 𝘅[2]
    #     d𝘅[3] = Κ[2] * 𝘅[2] - Κ[1]*𝘅[1]*𝘅[3]
    #     d𝘅[4] = Κ[3] * 𝘅[2] 
    # end

    function f!(𝘅,p,t)
        global specie, Κ
        d𝘅 = sum([𝛎[ℓ,i]*Jflux(𝘅,ℓ) for i in eachindex(specie), ℓ in eachindex(Κ)],dims=2)[:]
    end

    𝘅₀ = [2e-7; 0; 5e-7; 0]

    prob = ODEProblem(f!,𝘅₀,(T[1],T[end]))
    sol = solve(prob,RK4();dt= .5,adaptive=false)

    x = molecules ? hcat(sol.u...) * nₐ * V : hcat(sol.u...);

    mkpath(path)
    flname = path*"/"*model_nm;
    jldsave(flname, specie=specie, x=x, T=T)

    return (specie, x, T)
end

function Jflux(𝘅,ℓ)  # Concentration Flux
    return Κ[ℓ]*prod((^).(𝘅,Re[ℓ,:]))
end

function GibbsFreeEnergy(𝘅)    # Gibbs free energy
    g1 = [(-1)^ℓ * Jflux(𝘅,ℓ) for ℓ in eachindex(K)]
    g2 = [-(-1)^ℓ * log(Jflux(𝘅,ℓ)) for ℓ in eachindex(K)]
    g1[isinf.(g1) .|| g1 .== 0.0 ] .= 0.0;
    g2[isinf.(g2) .|| g2 .== 0.0 ] .= 0.0;
    g1g2 = sum(g1 * g2')
    d𝘅 = sum(g1g2[!isnan.(g1g2)]);
    return d𝘅
end