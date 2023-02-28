# using ProgressMeter
using StatsBase, Distributions

# function Gillespie(K, 𝛎, Re, S₀, T) # Gillespie
#     t = 0.0
#     iT = 1;
#     Sₜ = Vector{Int64}();
#     S = copy(S₀);
#     M = size(𝛎,1);

#     while t <= T[end]
#         if t >= T[iT]
#             for s in S
#                 push!(Sₜ,s)
#             end
#             iT += 1;
#         end

#         𝛂 = [K[m] * prod(α(S,Re,m)) for m in 1:M]
#         α₀ = sum(𝛂);
#         α₀ == 0.0 ? break : nothing

#         r = rand(Uniform(),1)
#         τ = log(1 / r[1]) / α₀;
#         t += τ
#         S += 𝛎[rand(Multinomial(1,vec(𝛂/α₀))) .!= 0,:];
#     end

#     if length(T)-iT > 0
#         for s in repeat(vec(S),length(T)-iT)
#             push!(Sₜ,s)
#         end
#     end

#     Sₜ = reshape(Sₜ,length(S₀),length(T)-1)';
#     return Sₜ
# end

# function Gillespie(K, 𝛎, Re, S₀, t, T) # Gillespie
#     N = size(Re,2);
#     Sₜ = copy(S₀);
#     S = copy(S₀);
#     M = size(𝛎,1);

#     α(𝓘,Re,m) = binomial.(𝓘,Re[m,:]') .* factorial.(Re[m,:]'); # Mass action law

#     while true
#         𝛂 = [K[m] * prod(α(S,Re,m)) for m in 1:M]
#         α₀ = sum(𝛂);
#         α₀ == 0.0 ? break : nothing

#         r = rand(Uniform(),2)

#         ir = rand(Multinomial(1,vec(𝛂/α₀))) .!= 0;
#         τ = log(1 / r[1]) / α₀;
        
#         if t+τ <= T[end]
#             t += τ
#             S += 𝛎[ir ,:];
#         else
#             break
#         end
#     end
#     return t, S
# end

function SSASolver(path, model_nm; saveprob=false, savestats=:eval)

    mkpath(path)

    model = "reactions/"*model_nm*".jl";
    include(model);

    # p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
    # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
    # # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
    # p₀ ./= sum(p₀);

    # # p₀ = ones(𝗻ₖ);              # Uniform distribution
    # # p₀ ./= sum(p₀); 
    # # p₀[end] = 1 - sum(p₀[1:end-1]);
    # cp(model,path*"/model.jl";force=true)
    
    # max_sim = 2500;
    # marg_labels = [];
    # marg = Vector{Any}(undef,length(T));
    # 𝔼 = zeros(length(𝗻ₖ),length(T));
    # # 𝕍ar = zeros(length(𝗻ₖ),length(T));
    # # Sk = zeros(length(𝗻ₖ),length(T));
    # 𝕊 = zeros(1,length(T));
    # # Si = zeros(1,length(T));
    # # Se = zeros(1,length(T));

    # uf = p₀[:];

    # if saveprob
    #     flname = path*"/"*model_nm*"_t"*string(0);
    #     jldsave(flname, p=uf, t=0)
    # end

    # println("Saving on "*path*".")
    
    # if (savestats != false)
    #     marg_labels, marg[1], 𝔼[:,1], 𝕍ar[:,1], Sk[:,1], 𝕊[1] = Statistics(uf,𝗻ₖ,specie);
    #     flname = path*"/"*model_nm*"_statistics_t"*string(0);
    #     jldsave(flname, specie=specie,marg_labels=marg_labels, 
    #         marg=marg[1], E=𝔼[:,1], Var=𝕍ar[:,1],Sk=Sk[:,1], S=𝕊[:,1], t=0, T=T)
    # end

    # pgres = Progress(length(T)-1; showspeed=true, desc="Solving with SSA...")

    # for iT in eachindex(T)[1:end-1]        
    #     pS = zeros(𝗻ₖ...);
    #     for _ in 1:max_sim

    #         𝒮 = [rand(ℰ),ℰ𝒜,rand(𝒜),ℬ]' .-1;
    #         S = Gillespie(K, 𝛎, Re, 𝒮, T[iT]);
    #         pS[(S .+ 1)...] += 1;
    #     end
    #     pS ./= sum(pS);

    #     if saveprob
    #         flname = path*"/"*model_nm*"_t"*string(iT);
    #         jldsave(flname, p=pS, t=T[iT+1])
    #     end

    #     if savestats == :eval
    #         marg_labels, marg[iT+1], 𝔼[:,iT+1], 𝕍ar[:,iT+1], Sk[:,iT+1], 𝕊[iT+1] = Statistics(pS,𝗻ₖ,specie);
    #     end

    #     if (savestats != false)
    #         flname = path*"/"*model_nm*"_statistics_t"*string(iT);
    #         jldsave(flname, specie=specie, marg_labels=marg_labels, 
    #         marg=marg[iT+1], E=𝔼[:,iT+1], Var=𝕍ar[:,iT+1],Sk=Sk[:,iT+1], S=𝕊[iT+1], t=T[iT+1], T=T)
    #     end

    #     ProgressMeter.next!(pgres)
    # end

    function Gillespie(K, 𝛎, Re, S₀, t, T) # Gillespie
        S = copy(S₀);
        M = size(𝛎,1);
        
        α(𝓘,Re,m,𝛎) = binomial.(𝓘,Re[m,:]') .* factorial.(Re[m,:]'); # Mass action law
        η(𝓘,Re,m,𝛎) = α(𝓘,Re,m,𝛎) .* ([0, 0, 0, 0] .<= (𝓘[:] + 𝛎[m,:]))' .* ((𝓘[:] + 𝛎[m,:]) .<= [(𝗻ₖ .- 1)...])';
    
        while t <= T[end]
            # 𝛂 = [K[m] * prod(η(S,Re,m,𝛎)) for m in 1:M]
            𝛂 = [K[m] * prod(α(S,Re,m,𝛎)) for m in 1:M]
            α₀ = sum(𝛂);
            α₀ == 0.0 ? break : nothing
    
            ir = rand(Multinomial(1,vec(𝛂/α₀))) .!= 0;
            τ = log(1 / rand()) / α₀;
            
            t += τ
            S += 𝛎[ir ,:];
        end
        return t, S
    end

    function BackwardsGillespie(K, 𝛎, Re, S₀, t, T) # Gillespie
        S = copy(S₀);
        M = size(𝛎,1);
        
        α(𝓘,Re,m,𝛎) = binomial.(𝓘 - 𝛎[m,:]',Re[m,:]') .* factorial.(Re[m,:]'); # Mass action law
        η(𝓘,Re,m,𝛎) = α(𝓘,Re,m,𝛎) .* ([0, 0, 0, 0] .<= (𝓘[:]))' .* ((𝓘[:]) .<= [(𝗻ₖ .- 1)...])';
    
        while t >= T[1]
            𝛂 = [K[m] * prod(η(S,Re,m,𝛎)) for m in 1:M]
            # 𝛂 = [K[m] * prod(α(S,Re,m,𝛎)) for m in 1:M]
            # println(S)
            α₀ = sum(𝛂);
            α₀ == 0.0 ? break : nothing
    
            ir = rand(Multinomial(1,vec(𝛂/α₀))) .!= 0;
            τ = log(1 / rand()) / α₀;
            
            t -= τ
            S -= 𝛎[ir ,:];
        end
        return t, S
    end
    
    Sent = zeros(length(T));
    E = zeros(length(T),length(𝗻ₖ));
    
    realizations = 2_000;
    𝒮 = (ℰ, ℰ𝒜, 𝒜, ℬ);
    R = hcat(rand.(map(x->x.-1 ,𝒮),realizations)...);
    TT = zeros(realizations);
    
    states = countmap(collect(eachrow(R)));
    p = collect(values(states))/realizations;
    Sent[1] = -sum(p .* log.(p));
    E[1,:] = sum(R, dims=1)/realizations; 

    pgres = Progress(length(T)-1; showspeed=true, desc="Solving the SSA...")
    
    for iT in eachindex(T)[2:end]
        # global p, states
        Threads.@threads for th in 1:realizations
            TT[th], R[th,:] = Gillespie(K,𝛎,Re, R[th,:]', TT[th], (T[iT-1],T[iT]))
        end
    
        # states = countmap(collect(eachrow(R)));
        # p = collect(values(states))/realizations;
        # Sent[iT] = -sum(p .* log.(p));
        # E[iT,:] = sum(R, dims=1)/realizations; 
        # # @printf "Evolution: %.0f%% \r" (t/length(T)*100)
        # ProgressMeter.next!(pgres)

    end

    TT = T[end] * ones(realizations);
    
    for iT in eachindex(T)[1:end-1]
        global p, states
        Threads.@threads for th in 1:realizations
            TT[th], R[th,:] = BackwardsGillespie(K,𝛎,Re, R[th,:]', TT[th], (T[end-iT],T[end-iT+1]))
        end

        states = countmap(collect(eachrow(R)));
        p = collect(values(states))/realizations;
        Sent[iT] = -sum(p .* log.(p));
        E[iT,:] = sum(R, dims=1)/realizations; 
        # @printf "Evolution: %.0f%% \r" (iT/length(T)*100)
        ProgressMeter.next!(pgres)

    end
    # flname = path*"/"*model_nm*"_statistics";
    # jldsave(flname, E=E, S=Sent, T=T)
    # return (marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊)
    return E, Sent
end
