using Pkg
Pkg.activate(".")
using ICRN

cod = "to58E_20230201"

path = "outputs/"*cod
model_nm = "MichaelisMenten"

saveStatistics(path*"/CME",model_nm);