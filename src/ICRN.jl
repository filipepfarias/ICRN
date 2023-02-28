module ICRN
    using MKL, MKLSparse, LinearAlgebra, SparseArrays

    include("CME.jl")

    # include("statistics.jl")

    # include("SSA.jl")

    # include("Deterministic.jl")

    # include("utils.jl")

    export operator, α, η
    # export CMEOperator, CMESolver, Statistics, Statistics!, saveStatistics, Entropy, 
    # dEntropy, Mean, Variance, Skewness, 
    # Marginals, Gillespie, SSASolver, DetSolver, Jflux, GibbsFreeEnergy
end