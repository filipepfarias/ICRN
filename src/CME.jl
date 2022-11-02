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
    W(𝓘,Re,m,𝛎) = reduce(kron,spdiagm.(eachcol(α(𝓘,Re,m))));
    H(𝓘,Re,m,𝛎) = reduce(kron,spdiagm.(eachcol(η(𝓘,Re,m,𝛎))));

    function CMEOperator(𝝼,Re,K,𝗻ₖ)
        𝓘 = hcat((:).(1,𝗻ₖ)...,);
        return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*W(𝓘,Re,m,𝝼) for m in eachindex(𝝼[:,1])]))
    end

    function CMEEntropy(X)
        return map(x -> -sum(filter(!isnan,x/sum(x) .* log.(x/sum(x)))),X)
    end

    function CMEMutualInformation(p,A,dt)
        
        ℐ = 0.0 .* p[1];
        for i in eachindex(p)[2:end]
            Pₖ₋₁ = repeat(p[i-1]', length(p[i-1]),1);
            Q = (I + A*dt) .* Pₖ₋₁;
            Pₖ = repeat(sum(Q,dims=2),1,length(p[i-1]));
            
            
            # Pₖ₋₁Pₖ = (I + A*dt) .* Pₖ₋₁;
            # Pₖ = sum(Pₖ₋₁Pₖ, dims=2);
            # Pₖ₋₁Pₖ[ Pₖ₋₁Pₖ .< 0 ] .= 0;
            Pₖ[ Pₖ .< 0 ] .= 0;
            Pₖ₋₁[ Pₖ₋₁ .< 0 ] .= 0;
            Q[ Q .< 0 ] .= 0;
            # Pₖ = repeat(Pₖ,1,length(Pₖ));
            ℐ[i] = sum( Q .* ( log.(Q) .- log.(Pₖ) .- log.(Pₖ₋₁) ) )
        end
        return ℐ
    end

    include("src/utils.jl")

    export CMEOperator, CMEEntropy, CMEMutualInformation
end