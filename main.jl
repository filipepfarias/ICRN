using Distributed
@everywhere using Pkg
@everywhere Pkg.activate(".")
using ICRN
using Random, Dates


path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")
mkpath(path)

model_nm = "MichaelisMenten"

sol_CME = CMESolver(path*"/CME", model_nm; saveprob=true, savestats=false)

# sol_SSA = SSASolver(path*"/SSA", model_nm; saveprob=false, savestats=:eval)

# sol_Det = DetSolver(path*"/Det", model_nm; molecules=true)