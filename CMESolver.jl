# using Pkg
# Pkg.activate(".")
# Pkg.instantiate();

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

T = 0.0:.2:50.0;
f(u,p,t) = A*u ;

begin
    u0 = p₀[:];
    p = u0;

    # flname = path*"/"*model_nm*"_t"*string(0);
    # jldsave(flname, p=p, t=0.0)
    # E = CMEMarginals(u0,𝗻ₖ,specie,flname,0.0)

    marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊 = CMEStatistics(u0,𝗻ₖ,specie);
    
    println("Saving on "*path*".")

    flname = path*"/"*model_nm*"_statistics_t"*string(0);
    jldsave(flname, specie=specie,
        marg_labels=marg_labels, 
        marg=marg, E=𝔼, Var=𝕍ar, R=ℝ, Sk=Sk, S=𝕊, t=0, T=T)

    pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

    for iT in eachindex(T)[1:end-1]
        global u0, flname, marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊
        prob = ODEProblem(f,u0, (T[iT],T[iT+1]));
        sol = solve(prob, RK4();dt= .5/15, saveat=T[iT+1],adaptive=false);
        u0 = sol.u[end]/sum(sol.u[end]);
        # u0 = sol.u[end]

        # flname = path*"/"*model_nm*"_t"*string(iT);
        # jldsave(flname, p=u0, t=T[iT+1])
        marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊 = CMEStatistics(u0,𝗻ₖ,specie)

        flname = path*"/"*model_nm*"_statistics_t"*string(iT);
        jldsave(flname, specie=specie,
        marg_labels=marg_labels, 
        marg=marg, E=𝔼, Var=𝕍ar, R=ℝ, Sk=Sk, S=𝕊, t=T[iT], T=T)

        ProgressMeter.next!(pgres)
    end
end

# Plotting
println("Saving plots...")
include("misc_plotting2.jl")
