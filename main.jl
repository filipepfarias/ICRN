using Pkg
Pkg.activate(".")
using ICRN
using Random, Dates


path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")
# path = "outputs/cnYSl_20230215"
mkpath(path)

model_nm = "MichaelisMenten"

sol_CME = CMESolver(path*"/CME", model_nm; saveprob=false, savestats=true)

sol_SSA = SSASolver(path*"/SSA", model_nm; saveprob=false, savestats=true)

# sol_Det = DetSolver(path*"/Det", model_nm; molecules=true)