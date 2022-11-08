using Pkg
Pkg.activate(".")
Pkg.instantiate();

using CME
using Random, Dates, FileIO, JLD2
using DifferentialEquations
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

T = 0.0:.5:50.0;
f(u,p,t) = A*u ;

begin
    u0 = p₀[:];

    p = u0;
    flname = path*"/MichaelisMenten_t"*string(0);
    jldsave(flname, p=p, t=0.0)
    
    println("Saving on "*path*".")
    pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

    for iT in eachindex(T)[1:end-1]
        global u0, flname
        prob = ODEProblem(f,u0, (T[iT],T[iT+1]));
        sol = solve(prob, RK4();dt=.5/2,saveat=T[iT+1],adaptive=false);
        sol.u[end][sol.u[end] .< 0] .= 0;
        # append!(p,[sol.u[end]]);
        u0 = sol.u[end]/sum(sol.u[end]);

        flname = path*"/MichaelisMenten_t"*string(iT);
        jldsave(flname, p=u0, t=T[iT+1])
        ProgressMeter.next!(pgres)
    end
end

specie = ["E", "EA", "A", "B"];
𝔼 = zeros(length(𝗻ₖ),size(T)...,);

pgres = Progress(length(T)-1; showspeed=true, desc="Computing statistics...")
for iT in eachindex(T)
    iT -= 1;
    local flname, p
    flname = path*"/MichaelisMenten_t"*string(iT);
    p = jldopen(flname)["p"];
    𝓅  = reshape(p,𝗻ₖ...,);
    𝓅ₙ = sum(𝓅);
    for i in 1:4, j in 1:4
        if j > i
            d = deleteat!(collect(1:4), [i j])
            mat = reshape(sum(𝓅,dims=d) ./ 𝓅ₙ ,𝗻ₖ[i],𝗻ₖ[j])'
            flsuffix = specie[i]*"_x_"*specie[j];
        elseif i == j
            ind = collect(1:4)
            mat = sum(𝓅,dims=deleteat!(ind,i))[:] ./𝓅ₙ ;
            flsuffix = specie[i];
        end
        if j >= i
            jldsave(flname*"_marg_"*flsuffix, p=mat, t=T[iT+1])
        end
    end

    𝔼[:,iT+1] = [
        sum(collect(0:(𝗻ₖ[i]-1)) .* sum(𝓅,dims=deleteat!(collect(1:length(𝗻ₖ)),i))[:] ./ 𝓅ₙ )
        for i in 1:length(𝗻ₖ)]
    
    ProgressMeter.next!(pgres)
end
flname = path*"/MichaelisMenten_mean";
jldsave(flname, E=𝔼)

# # Plotting
# println("Saving plots...")
# include("misc_plotting.jl")
