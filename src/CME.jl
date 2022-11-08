module CME
    using LinearAlgebra, SparseArrays

    function J(νi,n) # νi per reaction
        return νi > 0 ? sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] : sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
    end

    function 𝗝(ν,n)
        return reduce(kron,J.(ν,n))
    end

    α(𝓘,Re,m) = binomial.(𝓘,Re[m,:]')
    η(𝓘,Re,m,𝛎) = binomial.(𝓘,Re[m,:]') .* (𝓘 .<= (𝓘[end,:]' - 𝛎[m,:]')) .* (𝓘 .>= (𝓘[1,:]' - 𝛎[m,:]'));
    W(𝓘,Re,m,𝛎) = reduce(kron,Diagonal.(eachcol(α(𝓘,Re,m))));
    H(𝓘,Re,m,𝛎) = reduce(kron,Diagonal.(eachcol(η(𝓘,Re,m,𝛎))));

    function CMEOperator(𝝼,Re,K,𝗻ₖ)
        𝓘 = hcat((:).(1,𝗻ₖ)...,);
        return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*W(𝓘,Re,m,𝝼) for m in eachindex(𝝼[:,1])]));
    end

    function CMEEntropy(Xₖ)
        return -sum(filter(!isnan,Xₖ/sum(Xₖ) .* log.(Xₖ/sum(Xₖ))))
    end

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

    function CMERK(A,p,ΔT, N=3)
        δt = (ΔT[2] - ΔT[1])/N;
        k1 = k2 = k3 = k4 = similar(p);
        A *= δt; 
        for n in 1:N
            k1 = A*p;
            k2 = k1 + A*k1/2;
            k3 = k1 + A*k2/2;
            k4 = k1 + A*k3;

            p += (k1+2*k2+2*k3+k4)/6;
        end
        return p
    end

    export CMEOperator, CMEEntropy, CMEMutualInformation, CMERK
end