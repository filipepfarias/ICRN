using Distributed
addprocs()

@everywhere begin
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate();
    using CME
end

using Random, Dates, FileIO, JLD2
using DifferentialEquations: solve, ODEProblem, RK4
using ProgressMeter

path = "outputs/"*randstring(5)*"_"*Dates.format(now(),"yyyymmdd")
mkdir(path)

println("Building the CME operator...")
comp_time = @elapsed begin
    model_nm = "MichaelisMenten"
    model = "reactions/"*model_nm*".jl";
    include(model);
    A = CMEOperator(𝛎,Re,K,𝗻ₖ);   # CME Operator      
    cp(model,path*"/model.jl")
end
println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")

f(u,p,t) = A*u ;

uf = p₀[:];
p = uf;

marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊, Si, Se = CMEStatistics(uf,A,𝗻ₖ,specie);

println("Saving on "*path*".")

flname = path*"/"*model_nm*"_statistics_t"*string(0);
jldsave(flname, specie=specie,
    marg_labels=marg_labels, 
    marg=marg, E=𝔼, Var=𝕍ar, Sk=Sk, S=𝕊, Si=Si, Se=Se, t=0, T=T)

pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

@sync for iT in eachindex(T)[1:end-1]
    global pf, uf, flname, marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊, Si, Se, sol
for iT in eachindex(T)[1:end-1]
    global pf, uf, flname, marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊, Si, Se, sol
    prob = ODEProblem(f,uf, (T[iT],T[iT+1]));
    sol = solve(prob, RK4();dt= .5/20, saveat=T[iT+1],adaptive=false);
    uf = sol.u[end]/sum(sol.u[end]);
    # u0 = sol.u[end]
    # pf[:,iT+1] = uf;

    flname = path*"/"*model_nm*"_t"*string(iT);
    jldsave(flname, p=uf, t=T[iT+1])
    marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊, Si, Se = CMEStatistics(uf,A,𝗻ₖ,specie)

    flname = path*"/"*model_nm*"_statistics_t"*string(iT);
    @spawn jldsave(flname, specie=specie,
    marg_labels=marg_labels, 
    marg=marg, E=𝔼, Var=𝕍ar, Sk=Sk, S=𝕊, Si=Si, Se=Se, t=T[iT], T=T)

    ProgressMeter.next!(pgres)
end

# Plotting
# println("Saving plots...")
# include("misc_plotting.jl")
