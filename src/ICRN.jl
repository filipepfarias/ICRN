module ICRN
    # using  Distributed

    # @everywhere include("parallel_matmul.jl")

    include("CME.jl")

    include("statistics.jl")

    include("SSA.jl")

    include("Deterministic.jl")

    include("utils.jl")

    export CMEOperator, CMESolver, Statistics, Statistics!, Entropy, 
    dEntropy, Mean, Variance, Skewness, 
    Marginals, Gillespie, SSASolver, DetSolver, Jflux, GibbsFreeEnergy
end