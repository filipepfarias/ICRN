using Pkg
Pkg.activate(".")
using ICRN
using Random, Dates
# using Plots


path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")
# # path = "outputs/irxWK_20230223"
# mkpath(path)

model_nm = "Enzy_MichaelisMenten"
include("reactions/Enzy_MichaelisMenten.jl")

include("solve_CME.jl")

# sol_CME = CMESolver(path*"/CME", model_nm; saveprob=false, savestats=true)

# sol_SSA = SSASolver(path*"/SSA", model_nm; saveprob=false, savestats=true)

# sol_Det = DetSolver(path*"/Det", model_nm; molecules=true)

# plot(sol_CME[1]')

using Plots

