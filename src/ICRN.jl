module ICRN
    using LinearAlgebra, SparseArrays

    include("CME.jl")

    include("statistics.jl")

    # include("SSA.jl")

    # include("Deterministic.jl")

    # include("utils.jl")

    export operator, α, η, lotus, marginal, mean, entropy, KLdivergence
    # export CMEOperator, CMESolver, Statistics, Statistics!, saveStatistics, Entropy, 
    # dEntropy, Mean, Variance, Skewness, 
    # Marginals, Gillespie, SSASolver, DetSolver, Jflux, GibbsFreeEnergy
end