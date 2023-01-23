using ProgressMeter
using Distributed

@everywhere function Gillespie(K, 𝛎, Re, S₀, T) # Gillespie
    t = 0.0
    iT = 1;
    Sₜ = Vector{Int64}();
    S = copy(S₀);
    M = size(𝛎,1);

    while t <= T[end]
        if t >= T[iT]
            for s in S
                push!(Sₜ,s)
            end
            iT += 1;
        end

        𝛂 = [K[m] * prod(α(S,Re,m)) for m in 1:M]
        α₀ = sum(𝛂);
        α₀ == 0.0 ? break : nothing

        r = rand(Uniform(),1)
        τ = log(1 / r[1]) / α₀;
        t += τ
        S += 𝛎[rand(Multinomial(1,vec(𝛂/α₀))) .!= 0,:];
    end

    if length(T)-iT > 0
        for s in repeat(vec(S),length(T)-iT)
            push!(Sₜ,s)
        end
    end

    Sₜ = reshape(Sₜ,length(S₀),length(T)-1)';
    return Sₜ
end

function SSASolver(path, model_nm; saveprob=false, savestats=:eval)

    mkpath(path)

    @everywhere model = "reactions/"*model_nm*".jl";
    @everywhere include(model);

    p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
    p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
    # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
    p₀ ./= sum(p₀);

    # p₀ = ones(𝗻ₖ);              # Uniform distribution
    # p₀ ./= sum(p₀); 
    # p₀[end] = 1 - sum(p₀[1:end-1]);
    cp(model,path*"/model.jl";force=true)
    
    max_sim = 2500;
    marg_labels = [];
    marg = Vector{Any}(undef,length(T));
    𝔼 = zeros(length(𝗻ₖ),length(T));
    𝕍ar = zeros(length(𝗻ₖ),length(T));
    Sk = zeros(length(𝗻ₖ),length(T));
    𝕊 = zeros(1,length(T));
    Si = zeros(1,length(T));
    Se = zeros(1,length(T));

    uf = p₀[:];

    if saveprob
        flname = path*"/"*model_nm*"_t"*string(0);
        jldsave(flname, p=uf, t=0)
    end

    println("Saving on "*path*".")
    
    if (savestats != false)
        marg_labels, marg[1], 𝔼[:,1], 𝕍ar[:,1], Sk[:,1], 𝕊[1] = Statistics(uf,𝗻ₖ,specie);
        flname = path*"/"*model_nm*"_statistics_t"*string(0);
        jldsave(flname, specie=specie,marg_labels=marg_labels, 
            marg=marg[1], E=𝔼[:,1], Var=𝕍ar[:,1],Sk=Sk[:,1], S=𝕊[:,1], t=0, T=T)
    end

    pgres = Progress(length(T)-1; showspeed=true, desc="Solving with SSA...")

    for iT in eachindex(T)[1:end-1]        
        pS = zeros(𝗻ₖ...);
        for _ in 1:max_sim

            𝒮 = [rand(ℰ),ℰ𝒜,rand(𝒜),ℬ]' .-1;
            S = Gillespie(K, 𝛎, Re, 𝒮, T[iT]);
            pS[(S .+ 1)...] += 1;
        end
        pS ./= sum(pS);

        if saveprob
            flname = path*"/"*model_nm*"_t"*string(iT);
            jldsave(flname, p=pS, t=T[iT+1])
        end

        if savestats == :eval
            marg_labels, marg[iT+1], 𝔼[:,iT+1], 𝕍ar[:,iT+1], Sk[:,iT+1], 𝕊[iT+1] = Statistics(pS,𝗻ₖ,specie);
        end

        if (savestats != false)
            flname = path*"/"*model_nm*"_statistics_t"*string(iT);
            jldsave(flname, specie=specie, marg_labels=marg_labels, 
            marg=marg[iT+1], E=𝔼[:,iT+1], Var=𝕍ar[:,iT+1],Sk=Sk[:,iT+1], S=𝕊[iT+1], t=T[iT+1], T=T)
        end

        ProgressMeter.next!(pgres)
    end

    return (marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊)
end
