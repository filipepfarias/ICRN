module ICRN
    using MKL, MKLSparse, LinearAlgebra, SparseArrays

    include("CME.jl")

    include("statistics.jl")

    # include("SSA.jl")

    # include("Deterministic.jl")
    
    include("Macroscopic.jl")

    # include("utils.jl")

    export operatorα, operatorη, α, η, lotus, marginal, mean, entropy, 
    entropy_production, entropy_flow, KLdivergence, energy_input_rate,
    free_energy, d_entropy, sparselog, R, dxdt!
    # export CMEOperator, CMESolver, Statistics, Statistics!, saveStatistics, Entropy, 
    # dEntropy, Mean, Variance, Skewness, 
    # Marginals, Gillespie, SSASolver, DetSolver, Jflux, GibbsFreeEnergy
end