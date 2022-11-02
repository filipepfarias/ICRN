using CME, DifferentialEquations
using JLD2, FileIO
using Dates, Random

path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")
mkdir(path)

model = "reactions/MichaelisMenten.jl";
include(model);
cp(model,path*"/model.jl")

T = 0.0:.5:50.0;
f(u,p,t) = A*u ;

begin
    u0 = p₀[:];
    p = Vector{Vector{Float64}}([u0]);
    flname = path*"/MichaelisMenten_t"*string(0);
    jldsave(flname, p=p, t=0.0)
    for iT in eachindex(T)[1:end-1]
        global u0
        prob = ODEProblem(f,u0, (T[iT],T[iT+1]));
        sol = solve(prob, RK4(),saveat=T[iT+1]);
        sol.u[end][sol.u[end] .< 0] .= 0;
        append!(p,[sol.u[end]]);
        u0 = sol.u[end]/sum(sol.u[end]);
        flname = path*"/MichaelisMenten_t"*string(iT);
        jldsave(flname, p=u0, t=T[iT+1])
    end
end
