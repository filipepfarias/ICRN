using LinearAlgebra, SparseArrays
using DifferentialEquations: solve, ODEProblem, RK4
using FileIO, JLD2
using ProgressMeter

function J(νi,n) # νi per reaction
    return νi > 0 ? sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] : sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
end

function 𝗝(ν,n)
    return kron(reverse(J.(ν,n))...)
end

vecoper(f,x,y) = map((x,y) -> f.(x,y),x,y)

α(𝓘,Re,m) = vecoper(binomial,𝓘,Re[m,:]) .* factorial.(Re[m,:]);
function η(𝓘,Re,m,𝛎) 
    ν1 = vecoper(-,getindex.(𝓘,length.(𝓘)),𝛎[m,:]);
    ν1 = vecoper(<=,𝓘,ν1);

    ν2 = vecoper(-,getindex.(𝓘,(ones(Int,length(𝓘))...,)),𝛎[m,:]);
    ν2 = vecoper(>=,𝓘,ν2);

    res = vecoper(*,ν1,ν2);
    res = vecoper(*,α(𝓘,Re,m),res);
    return res
end
W(𝓘,Re,m,𝛎) = reduce(kron,reverse(Diagonal.(α(𝓘,Re,m))));
H(𝓘,Re,m,𝛎) = reduce(kron,reverse(Diagonal.(η(𝓘,Re,m,𝛎))));

function CMEOperator(𝝼,Re,K,𝗻ₖ)
    𝓘 = [collect.((:).(1,𝗻ₖ))...,];
    return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*H(𝓘,Re,m,𝝼) for m in eachindex(𝝼[:,1])]));
end

function CMESolver(path, model_nm; saveprob=false, savestats=:eval)

    mkpath(path)
    println("Building the CME operator...")
    comp_time = @elapsed begin
        model = "reactions/"*model_nm*".jl";
        include(model);

        p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
        p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
        # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
        p₀ ./= sum(p₀);

        # p₀ = ones(𝗻ₖ);              # Uniform distribution
        # p₀ ./= sum(p₀); 
        # p₀[end] = 1 - sum(p₀[1:end-1]);
        A = CMEOperator(𝛎,Re,K,𝗻ₖ);   # CME Operator      
        cp(model,path*"/model.jl")
    end
    println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")
    
    marg_labels = [];
    marg = Vector{Any}(undef,length(T));
    𝔼 = zeros(length(𝗻ₖ),length(T));
    𝕍ar = zeros(length(𝗻ₖ),length(T));
    Sk = zeros(length(𝗻ₖ),length(T));
    𝕊 = zeros(1,length(T));
    Si = zeros(1,length(T));
    Se = zeros(1,length(T));


    f(u,p,t) = A*u ;

    uf = p₀[:];

    if saveprob
        flname = path*"/"*model_nm*"_t"*string(0);
        jldsave(flname, p=uf, t=0)
    end

    println("Saving on "*path*".")
    
    if (savestats != false)
        marg_labels, marg[1], 𝔼[:,1], 𝕍ar[:,1], Sk[:,1], 𝕊[1], Si[1], Se[1] = Statistics(uf,A,𝗻ₖ,specie);
        flname = path*"/"*model_nm*"_statistics_t"*string(0);
        jldsave(flname, specie=specie,marg_labels=marg_labels, 
            marg=marg[1], E=𝔼[:,1], Var=𝕍ar[:,1],Sk=Sk[:,1], S=𝕊[:,1], Si=Si[:,1], Se=Se[:,1], t=0, T=T)
    end

    pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

    for iT in eachindex(T)[1:end-1]        
        prob = ODEProblem(f,uf, (T[iT],T[iT+1]));
        sol = solve(prob, RK4();dt= .5/20, saveat=T[iT+1],adaptive=false);
        uf = sol.u[end]/sum(sol.u[end]);

        if saveprob
            flname = path*"/"*model_nm*"_t"*string(iT);
            jldsave(flname, p=uf, t=T[iT+1])
        end

        if savestats == :eval
            marg_labels, marg[iT+1], 𝔼[:,iT+1], 𝕍ar[:,iT+1], Sk[:,iT+1], 𝕊[iT+1], Si[iT+1], Se[iT+1] = Statistics(uf,A,𝗻ₖ,specie);
        end

        if (savestats != false)
            flname = path*"/"*model_nm*"_statistics_t"*string(iT);
            jldsave(flname, specie=specie, marg_labels=marg_labels, 
            marg=marg[iT+1], E=𝔼[:,iT+1], Var=𝕍ar[:,iT+1],Sk=Sk[:,iT+1], S=𝕊[iT+1], Si=Si[iT+1], Se=Se[iT+1], t=T[iT+1], T=T)
        end

        ProgressMeter.next!(pgres)
    end

    return (marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊, Si, Se)
end