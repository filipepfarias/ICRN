module CME
    using LinearAlgebra, SparseArrays

    function J(νi,n) # νi per reaction
        return νi > 0 ? sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] : sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
    end

    function 𝗝(ν,n)
        return kron(reverse(J.(ν,n))...)
    end

    α(𝓘,Re,m) = binomial.(𝓘,Re[m,:]') .* factorial.(Re[m,:]');
    η(𝓘,Re,m,𝛎) = α(𝓘,Re,m) .* (𝓘 .<= (𝓘[end,:]' - 𝛎[m,:]')) .* (𝓘 .>= (𝓘[1,:]' - 𝛎[m,:]'));
    W(𝓘,Re,m,𝛎) = reduce(kron,reverse(Diagonal.(eachcol(α(𝓘,Re,m)))));
    H(𝓘,Re,m,𝛎) = reduce(kron,reverse(Diagonal.(eachcol(η(𝓘,Re,m,𝛎)))));

    function CMEOperator(𝝼,Re,K,𝗻ₖ)
        𝓘 = hcat((:).(1,𝗻ₖ)...,);
        return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*W(𝓘,Re,m,𝝼) for m in eachindex(𝝼[:,1])]));
    end

    # function CMEEntropy(p,A)
    #     p = p .* log.(p)
    #     p[isnan.(p)] .= 0.0;
    #     return (-sum(p),-sum(A*p))
    # end

    function CMEMutualInformation(Xₖ₋₁,Xₖ,A,dt)
        
        ℐ = 0.0;

        for i in eachindex(Xₖ), j in eachindex(Xₖ₋₁)
            expA = (I[i,j]+A[i,j]*dt);
            logXₖₖ₋₁ = log(expA*Xₖ₋₁[j]);
            logXₖ₋₁  = log(Xₖ₋₁[j]);
            logXₖ   = log(Xₖ[i]);
            
            if !any(isinf.([logXₖₖ₋₁, logXₖ₋₁, logXₖ]))
                ℐ += expA*Xₖ₋₁[j] * (logXₖₖ₋₁ - logXₖ₋₁ - logXₖ)
            end
        end

        return ℐ
    end

    function CMEStatistics(p,A,𝗻ₖ,specie)
        𝓅  = reshape(p,𝗻ₖ...,);
        𝓅ₙ = sum(𝓅);

        𝕊 = dropzeros(p) .* log.(dropzeros(p))
        𝕊 = -sum(𝕊);

        Q = A - spdiagm(diag(A));
        logA = copy(Q);
        nzlogA = nonzeros(logA); nzlogA .= log.(nonzeros(logA));

        logAp = Q*dropzeros(spdiagm(p));
        nzlogAp = nonzeros(logAp); nzlogAp .= log.(nonzeros(logAp));

        J = Q*dropzeros(spdiagm(p)) - (Q*dropzeros(spdiagm(p)))';
        X = logAp - logAp';

        Se = .5 * sum(J .* (logA - logA')); # h_ex
        Si = .5 * sum( J .* X ); # e_p

        # Si = .5* sum( (A*spdiagm(p) - (A*spdiagm(p))') .* (log.(Q * spdiagm(p)) .* ((Q * spdiagm(p)) .!= 0) - log.(Q * spdiagm(p))' .* ((Q * spdiagm(p)) .!= 0)' ) )
        # Se = .5* sum( (A*spdiagm(p) - (A*spdiagm(p))') .* (log.(Q) .* (Q .!= 0) - log.(Q') .* (Q' .!= 0) ) ) 

        𝔼 = zeros(length(𝗻ₖ),1);
        𝕍ar = zeros(length(𝗻ₖ),1);
        ℝ = zeros(length(𝗻ₖ),1,length(𝗻ₖ));
        Sk = zeros(length(𝗻ₖ),1);

        idn = Int(length(𝗻ₖ)*(length(𝗻ₖ)+1)/2);
        marg = Vector{Any}(undef,idn);
        idm = 1;
        marg_labels = [];

        for i in eachindex(𝗻ₖ), j in eachindex(𝗻ₖ)
            if j > i
                ℕxℕ = ((1:𝗻ₖ[i]) .- 1)*((1:𝗻ₖ[i]) .- 1)';
                d = deleteat!(collect(eachindex(𝗻ₖ)), [i j])
                marg[idm] = reshape(sum(𝓅,dims=d)./𝓅ₙ ,𝗻ₖ[i],𝗻ₖ[j])

                ℝ[i,1,j] = sum( ℕxℕ .* marg[idm] );

                flsuffix = specie[i]*"_x_"*specie[j];
                append!(marg_labels,[flsuffix]);
                idm += 1;
            elseif i == j
                ℕ = (1:𝗻ₖ[i]) .- 1;
                ind = collect(1:length(𝗻ₖ))
                marg[idm] = sum(𝓅,dims=deleteat!(ind,i))[:] ./𝓅ₙ ;
                𝔼[i,1] = sum( ℕ .* marg[idm] );
                𝕍ar[i,1] = sum( (ℕ.-𝔼[i,1]).^2 .* marg[idm] );
                Sk[i,1] = 𝕍ar[i,1] == 0.0 ? 0.0 : sum(((ℕ.-𝔼[i,1])./√𝕍ar[i,1]).^3 .* marg[idm] );

                flsuffix = specie[i];
                append!(marg_labels,[flsuffix]);
                idm += 1;
            end
        end

        return marg_labels, marg, 𝔼, 𝕍ar, ℝ, Sk, 𝕊, Si, Se
    end

    function Gillespie(K, 𝛎, S₀, T) # Gillespie
        t = 0
        t_vec = [0.0]
        Sₜ = S₀;
        S = S₀;
        # δ = [-2 -1 1 0; 0 -1 0 +1]';
        d = size(𝛎,2);
    
        while t <= T
            # α = k .* [Sₜ[1]*(Sₜ[1]-1) Sₜ[1]*Sₜ[2] 1 1]
            𝛂 = [K[m] * α(Sₜ,Re,m) for m in 1:d]
            α₀ = sum(𝛂);
            # αᵣ = cumsum(α,dims=1)/α₀;
        
            r = rand(Uniform(),1)
            τ = log(1 / r[1]) / α₀;
            t += τ
            append!(t_vec,t)
    
            Sₜ += 𝛎[rand(Multinomial(1,vec(𝛂/α₀))) .!= 0,:];
            S = cat(S,Sₜ,dims=1)
        end
        return t_vec,S
    end

    export CMEOperator, CMEStatistics, Gillespie
end