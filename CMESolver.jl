using Pkg
Pkg.activate(".")
Pkg.instantiate();

using CME
using Random, Dates, FileIO, JLD2
using DifferentialEquations: solve, ODEProblem, RK4
using ProgressMeter

path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")
mkdir(path)

println("Building the CME operator...")
comp_time = @elapsed begin
    model = "reactions/MichaelisMenten.jl";
    include(model);
    cp(model,path*"/model.jl")
end
println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")

T = 0.0:1:100.0;
f(u,p,t) = A*u ;

begin
    u0 = p₀[:];
    p = u0;
    𝔼 = zeros(length(𝗻ₖ),size(T)...,);
    ℍ = zeros(1,size(T)...,);

    flname = path*"/MichaelisMenten_t"*string(0);
    # jldsave(flname, p=p, t=0.0)
    E = CMEMarginals(u0,𝗻ₖ,specie,flname,0.0)
    
    𝔼[:,1] = E;
    ℍ[1] = CMEEntropy(u0);
    
    println("Saving on "*path*".")
    pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

    for iT in eachindex(T)[1:end-1]
        global u0, flname, 𝔼, E, ℍ
        prob = ODEProblem(f,u0, (T[iT],T[iT+1]));
        sol = solve(prob, RK4();dt=.5/4,saveat=T[iT+1],adaptive=false);
        # sol.u[end][sol.u[end] .< 0] .= 0;
        # append!(p,[sol.u[end]]);
        u0 = sol.u[end]/sum(sol.u[end]);
        # u0 = sol.u[end]

        flname = path*"/MichaelisMenten_t"*string(iT);
        # jldsave(flname, p=u0, t=T[iT+1])
        E = CMEMarginals(u0,𝗻ₖ,specie,flname,T[iT])

        𝔼[:,iT+1] = E;
        # ℍ[iT+1] = CMEEntropy(u0);
        ProgressMeter.next!(pgres)
    end
    flname = path*"/MichaelisMenten_mean";
    jldsave(flname, E=𝔼)

    flname = path*"/MichaelisMenten_entropy";
    jldsave(flname, H=ℍ)
end


# pgres = Progress(length(T)-1; showspeed=true, desc="Computing statistics...")
# for iT in eachindex(T)
#     iT -= 1;
#     local flname, p, mat
#     flname = path*"/MichaelisMenten_t"*string(iT);
#     p = jldopen(flname)["p"];
#     rm(flname)
#     CMEMarginals(p,𝗻ₖ,specie,flname,T[iT+1])

    
    
#     ProgressMeter.next!(pgres)
# end
# flname = path*"/MichaelisMenten_mean";
# jldsave(flname, E=𝔼)

# Plotting
println("Saving plots...")
include("misc_plotting.jl")
