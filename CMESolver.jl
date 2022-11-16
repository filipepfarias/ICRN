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
    model_nm = "2ABreaction"
    model = "reactions/"*model_nm*".jl";
    include(model);
    cp(model,path*"/model.jl")
end
println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")

T = 0.0:1:200.0;
f(u,p,t) = A*u ;

begin
    u0 = p₀[:];
    p = u0;
    𝔼 = zeros(length(𝗻ₖ),size(T)...,);
    𝕊 = zeros(1,size(T)...,);
    d𝕊dt = zeros(1,size(T)...,);

    flname = path*"/"*model_nm*"_t"*string(0);
    jldsave(flname, p=p, t=0.0)
    E = CMEMarginals(u0,𝗻ₖ,specie,flname,0.0)
    
    𝔼[:,1] = E;
    𝕊[1],d𝕊dt[1] = CMEEntropy(u0,A);
    
    println("Saving on "*path*".")
    pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

    for iT in eachindex(T)[1:end-1]
        global u0, flname, 𝔼, E, 𝕊, d𝕊dt
        prob = ODEProblem(f,u0, (T[iT],T[iT+1]));
        sol = solve(prob, RK4();dt= .5/8, saveat=T[iT+1],adaptive=false);
        u0 = sol.u[end]/sum(sol.u[end]);
        # u0 = sol.u[end]

        flname = path*"/"*model_nm*"_t"*string(iT);
        jldsave(flname, p=u0, t=T[iT+1])
        E = CMEMarginals(u0,𝗻ₖ,specie,flname,T[iT])

        𝔼[:,iT+1] = E;
        𝕊[iT+1], d𝕊dt[iT+1] = CMEEntropy(u0,A);
        ProgressMeter.next!(pgres)
    end
    flname = path*"/"*model_nm*"_mean";
    jldsave(flname, E=𝔼)

    flname = path*"/"*model_nm*"_entropy";
    jldsave(flname, S=𝕊, dSdt=d𝕊dt)
end

# Plotting
println("Saving plots...")
include("misc_plotting.jl")
