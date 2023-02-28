using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using MKL, MKLSparse, SparseArrays, LinearAlgebra
using DifferentialEquations

# mkpath(path)
println("Building the CME operator...")
comp_time = @elapsed begin
    model = "reactions/"*model_nm*".jl";
    include(model);

    p₀ = spzeros(prod(𝗻ₖ));                # Initial condition for Section 7.3
    cart2lin = LinearIndices((:).(1,𝗻ₖ))[ℰ, ℰ𝒜, 𝒜, ℬ][:];
    p₀[cart2lin] .= 1.0;                  # Uniform distribution
    
    nzp₀ = nonzeros(p₀); nzp₀ .= 1/sum(nonzeros(p₀)); # Renormalizing distribution

    A = operator(𝛎,Re,K,𝗻ₖ);    # CME Operator 
    # cp(model,path*"/"*model_nm*".jl"; force=true)
end
println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")

# marg_labels = [];
# marg = Vector{Any}(undef,length(T));
𝔼 = zeros(length(𝗻ₖ),length(T));
# 𝕍ar = zeros(length(𝗻ₖ),length(T));
# Sk = zeros(length(𝗻ₖ),length(T));
𝕊 = zeros(1,length(T));
# Si = zeros(1,length(T));
# Se = zeros(1,length(T));

function forwardKolmogorov(dpf,pf,params,t) 
    if t == params;
        pf = pf/sum(pf);
    end
    mul!(dpf,A,pf)
    # dpf = A * pf; 
    return dpf
end

function backwardKolmogorov(dpb,pb,params,t) 
    if t == params;
        pb = pb/sum(pb);
    end
    mul!(dpb,pb,A) 
    # dpb = pb * A; 
    return dpb
end

println("Saving on "*path*".")
println("Solving the CME.")

pf = p₀[:]; # Forward probability initial distribution
    
fw_sol = solve(ODEProblem(forwardKolmogorov,pf,nothing,T[2]),
            RK4(); 
            dt= .5/10, 
            tspan = (T[1],T[end]),
            saveat=T[2],
            adaptive=false, progress=true, progress_steps = 10);

pb = sparse(fw_sol.u[end]');

bw_sol = solve(ODEProblem(backwardKolmogorov,ub,nothing,T[2]),
            RK4(); 
            dt= .5/10, 
            tspan = (T[1],T[end]),
            saveat=T[2],
            adaptive=false, progress=true, progress_steps = 10);
