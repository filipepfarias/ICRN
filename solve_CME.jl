using ProgressMeter
using FileIO, JLD2

using MKLSparse, SparseArrays, LinearAlgebra
using DifferentialEquations

mkpath(path)
println("Building the CME operator...")
comp_time = @elapsed begin
    model = "reactions/"*model_nm*".jl";
    include(model);

    p₀ = spzeros(prod(𝗻ₖ));                # Initial condition for Section 7.3
    # p₀ = zeros(prod(𝗻ₖ));                # Initial condition for Section 7.3
    cart2lin = LinearIndices((:).(1,𝗻ₖ))[ℰ, ℰ𝒜, 𝒜, ℬ][:];
    p₀[cart2lin] .= 1.0;                  # Uniform distribution
    
    nzp₀ = nonzeros(p₀); nzp₀ .= 1/sum(nonzeros(p₀)); # Renormalizing distribution
    # p₀ /= sum(p₀);

    A = operator(𝛎,Re,K,𝗻ₖ);    # CME Operator 
    At = sparse(A');
    cp(model,path*"/"*model_nm*".jl"; force=true)
end
println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")


function f(u,p,t) 
    return A*u;
end

function g(u,p,t) 
    return At*u;
end

println("Saving on "*path*".")

pf = copy(p₀); # Forward probability
pf_log = Vector{Any}();
push!(pf_log,pf);

pgres = Progress(length(T)-1; showspeed=true, desc="Solving the Forward CME...")
for iT in eachindex(T)[1:end-1] 
    global pf, pf_log    
    prob = ODEProblem(f,pf, (T[iT],T[iT+1]));
    fw_sol = solve(prob, RK4();dt= .5/20, saveat=T[iT+1],adaptive=false);
    pf = fw_sol.u[end]/sum(fw_sol.u[end]);
    push!(pf_log,pf);
    ProgressMeter.next!(pgres)
end

pb = copy(pf);
pb_log = Vector{Any}();
push!(pb_log,pb);

pgres = Progress(length(T)-1; showspeed=true, desc="Solving the Backward CME...")
for iT in eachindex(T)[1:end-1] 
    global pb, pb_log     
    prob = ODEProblem(g,pb, (T[iT],T[iT+1]));
    bw_sol = solve(prob, RK4();dt= .5/20, saveat=T[iT+1],adaptive=false);
    pb = bw_sol.u[end]/sum(bw_sol.u[end]);
    push!(pb_log,pb);
    ProgressMeter.next!(pgres)
end

jldsave(path*"/CME/"*model_nm, pf_log=pf_log, pb_log=pb_log);