using Pkg
Pkg.activate(".")
using ICRN

cod = "qyfWO_20230130"

path = "outputs/"*cod
model_nm = "MichaelisMenten"

saveStatistics(path*"/CME",model_nm);