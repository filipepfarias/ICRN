"""
J(νi,n)

Constructs the shift matrix for the νi-th reaction and the n-th species.
Thanks to Minho.
"""
function J(νi,n)
    if νi > 0
        return sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] 
    else 
        return sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
    end
    # if νi > 0
    #     return GBMatrix(collect(1:(n+νi)),collect(1:(n+νi)),true, n+νi, n+νi)[1:end-νi,νi+1:end] 
    # else 
    #     return GBMatrix(collect(1:(n+νi)),collect(1:(n+νi)),true, n-νi, n-νi)[1-νi:end,1:(end+νi)]
    # end
end

"""
𝗝(ν,n)

Constructs the kronecker product of shift matrices constructed by J(νi,n), 
using the `reverse` to force an lexicographic ordering of the product 
(the Kronecker product is not commutative).
"""
function 𝗝(ν,n)
    return kron(reverse(J.(ν,n))...)
end

"""
(..)(f,(x,y))

Perform `f` broadcasted between uneven vectors `x` and `y`.
"""
(..)(f,(x,y)) = map((x,y) -> f.(x,y),x,y)

"""
α(𝓘,Re,m)

Evaluate the propensity for 𝓘 species. Re is the stoichiometric coefficients 
of the reactants in the m-th reaction.
"""
α(𝓘,Re) = ..(binomial,(𝓘,Re)) .* factorial.(Re);

"""
η(𝓘,Re,m,𝛎)

Evaluate the propensity for 𝓘 species in the m-th reaction, but giving 0 
to the reaction that goes outside the state-space. Re is the stoichiometric 
coefficients of the reactants in the m-th reaction.
"""
function η(𝓘,Re,𝛎) 
    ν1 = ..(<=,(𝓘,..(-,(getindex.(𝓘,length.(𝓘)),𝛎))));
    ν2 = ..(>=,(𝓘,..(-,(getindex.(𝓘,(ones(Int,length(𝓘))...,)),𝛎))));

    res = ..(*,(ν1,ν2));
    res = ..(*,(α(𝓘,Re),res));
    return res
end

"""
diagα(𝓘,Re,m,𝛎)

Creates the diagonal matrix with the propensities for the m-th reaction[^1] 
using α(𝓘,Re).

[^1]: See Eq. 7 in [An adaptive solution to the chemical master equation using 
quantized tensor trains with sliding windows, Dinh and Sidje](https://doi.org/10.1088/1478-3975/aba1d2).
"""
diagα(𝓘,Re,m,𝛎) = reduce(kron,sparse.(reverse(Diagonal.(α(𝓘,Re[m,:])))));

"""
diagη(𝓘,Re,m,𝛎)

Creates the diagonal matrix with the propensities for the m-th reaction[^1] 
using η(𝓘,Re).

[^1]: See Eq. 7 in [An adaptive solution to the chemical master equation using 
quantized tensor trains with sliding windows, Dinh and Sidje](https://doi.org/10.1088/1478-3975/aba1d2).
"""
diagη(𝓘,Re,m,𝛎) = reduce(kron,sparse.(reverse(Diagonal.(η(𝓘,Re[m,:],𝛎[m,:])))));

"""
operatorη(𝛎,Re,K,𝗻ₖ)

Builds the trasition rate matrix for the chemical master equation[^1].

[^1]: See Eq. 7 in [An adaptive solution to the chemical master equation using 
quantized tensor trains with sliding windows, Dinh and Sidje](https://doi.org/10.1088/1478-3975/aba1d2).
"""
function operatorη(𝛎,Re,K,𝗻ₖ)
    𝓘 = [collect.((:).(0,𝗻ₖ .-1))...,];
    return operator(𝛎,Re,K,𝗻ₖ,𝓘);
end

"""
operatorη(𝛎,Re,K,𝗻ₖ,𝓘)

Builds the trasition rate matrix for the chemical master equation[^1].

[^1]: See Eq. 7 in [An adaptive solution to the chemical master equation using 
quantized tensor trains with sliding windows, Dinh and Sidje](https://doi.org/10.1088/1478-3975/aba1d2).
"""
function operatorη(𝛎,Re,K,𝗻ₖ,𝓘)
    return sum([(𝗝(𝛎[m,:],𝗻ₖ) - I)*K[m]*diagη(𝓘,Re,m,𝛎) for m in eachindex(𝛎[:,1])]);
end

"""
operatorα(𝛎,Re,K,𝗻ₖ)

Builds the trasition rate matrix for the chemical master equation[^1].

[^1]: See Eq. 7 in [An adaptive solution to the chemical master equation using 
quantized tensor trains with sliding windows, Dinh and Sidje](https://doi.org/10.1088/1478-3975/aba1d2).
"""
function operatorα(𝛎,Re,K,𝗻ₖ)
    𝓘 = [collect.((:).(0,𝗻ₖ .-1))...,];
    return operatorα(𝛎,Re,K,𝗻ₖ,𝓘);
end

"""
operatorα(𝛎,Re,K,𝗻ₖ,𝓘)

Builds the trasition rate matrix for the chemical master equation[^1].

[^1]: See Eq. 7 in [An adaptive solution to the chemical master equation using 
quantized tensor trains with sliding windows, Dinh and Sidje](https://doi.org/10.1088/1478-3975/aba1d2).
"""
function operatorα(𝛎,Re,K,𝗻ₖ,𝓘)
    return sum([(𝗝(𝛎[m,:],𝗻ₖ) - I)*K[m]*diagα(𝓘,Re,m,𝛎) for m in eachindex(𝛎[:,1])]);
end

################################## OLD

# using LinearAlgebra, SparseArrays, MKLSparse
# using FileIO, JLD2
# using ProgressMeter

# function J(νi,n) # νi per reaction
#     return νi > 0 ? sparse(I,n+νi,n+νi)[1:end-νi,νi+1:end] : sparse(I,n-νi,n-νi)[1-νi:end,1:(end+νi)]
# end

# function 𝗝(ν,n)
#     return kron(reverse(J.(ν,n))...)
# end

# vecoper(f,x,y) = map((x,y) -> f.(x,y),x,y)

# α(𝓘,Re,m) = vecoper(binomial,𝓘,Re[m,:]) .* factorial.(Re[m,:]);
# function η(𝓘,Re,m,𝛎) 
#     ν1 = vecoper(-,getindex.(𝓘,length.(𝓘)),𝛎[m,:]);
#     ν1 = vecoper(<=,𝓘,ν1);

#     ν2 = vecoper(-,getindex.(𝓘,(ones(Int,length(𝓘))...,)),𝛎[m,:]);
#     ν2 = vecoper(>=,𝓘,ν2);

#     res = vecoper(*,ν1,ν2);
#     res = vecoper(*,α(𝓘,Re,m),res);
#     return res
# end
# W(𝓘,Re,m,𝛎) = reduce(kron,reverse(Diagonal.(α(𝓘,Re,m))));
# H(𝓘,Re,m,𝛎) = reduce(kron,reverse(Diagonal.(η(𝓘,Re,m,𝛎))));

# function CMEOperator(𝝼,Re,K,𝗻ₖ)
#     𝓘 = [collect.((:).(0,𝗻ₖ .-1))...,];
#     return (sum([(𝗝(𝝼[m,:],𝗻ₖ) - I)*K[m]*W(𝓘,Re,m,𝝼) for m in eachindex(𝝼[:,1])]));
# end

# function CMESolver(path, model_nm; saveprob=false, savestats=:eval)

#     mkpath(path)
#     println("Building the CME operator...")
#     comp_time = @elapsed begin
#         model = "reactions/"*model_nm*".jl";
#         include(model);

#         p₀ = zeros(𝗻ₖ);                # Initial condition for Section 7.3
#         p₀[ℰ, ℰ𝒜, 𝒜, ℬ] .= 1.0;
#         # p₀[ℰ, ℰ𝒜, 𝒜, ℬ] = 1.0;
#         p₀ ./= sum(p₀);

#         # p₀ = ones(𝗻ₖ);              # Uniform distribution
#         # p₀ ./= sum(p₀); 
#         # p₀[end] = 1 - sum(p₀[1:end-1]);
#         precompile(CMEOperator,(Matrix{Int64},Matrix{Int64},Vector{Float64},NTuple{4,Int64}));
#         A = CMEOperator(𝛎,Re,K,𝗻ₖ);  # CME Operator 
#         cp(model,path*"/model.jl")
#     end
#     println("Computation time for the assemble of the operator: "*string(comp_time)*"s.")
    
#     # marg_labels = [];
#     # marg = Vector{Any}(undef,length(T));
#     𝔼 = zeros(length(𝗻ₖ),length(T));
#     # 𝕍ar = zeros(length(𝗻ₖ),length(T));
#     # Sk = zeros(length(𝗻ₖ),length(T));
#     𝕊 = zeros(1,length(T));
#     # Si = zeros(1,length(T));
#     # Se = zeros(1,length(T));

#     # function f(u,p,t) 
#     #     nt = BLAS.get_num_threads()
#     #     BLAS.set_num_threads(1)
#     #     F = A*u 
#     #     BLAS.set_num_threads(nt)
#     #     return F
#     # end
#     function f(u,p,t) 
#         return A*u 
#     end

#     function g(u,p,t) 
#         return u'*A 
#     end

#     uf = p₀[:];

#     # if saveprob
#     #     flname = path*"/"*model_nm*"_t"*string(0);
#     #     jldsave(flname, p=uf, t=0; compress = true)
#     # end

#     println("Saving on "*path*".")
    
#     # if (savestats != false)
#     #     marg_labels, marg[1], 𝔼[:,1], 𝕍ar[:,1], Sk[:,1], 𝕊[1], Si[1], Se[1] = Statistics(uf,A,𝗻ₖ,specie);
#     #     flname = path*"/"*model_nm*"_statistics_t"*string(0);
#     #     jldsave(flname, specie=specie,marg_labels=marg_labels, 
#     #         marg=marg[1], E=𝔼[:,1], Var=𝕍ar[:,1],Sk=Sk[:,1], S=𝕊[:,1], Si=Si[:,1], Se=Se[:,1], t=0, T=T)
#     # end

#     # 𝕊[1] = Entropy(uf);
#     # E = [sum(reshape(uf,𝗻ₖ...), dims=deleteat!(collect(1:length(𝗻ₖ)),i))[:] for i in 1:length(𝗻ₖ)];
#     # E = sum.(map(.*,E,collect.((:).(0,𝗻ₖ .- 1)))); 
#     # 𝔼[:,1]=E;

#     pgres = Progress(length(T)-1; showspeed=true, desc="Solving the CME...")

#     for iT in eachindex(T)[1:end-1]        
#         prob = ODEProblem(f,uf, (T[iT],T[iT+1]));
#         sol = solve(prob, RK4();dt= .5/20, saveat=T[iT+1],adaptive=false);
#         uf = sol.u[end]/sum(sol.u[end]);

#         # if saveprob
#         #     flname = path*"/"*model_nm*"_t"*string(iT);
#         #     jldsave(flname, p=uf, t=T[iT+1]; compress = true)
#         # end

#         # if savestats 
#         #     marg_labels, marg[iT+1], 𝔼[:,iT+1], 𝕍ar[:,iT+1], Sk[:,iT+1], 𝕊[iT+1], Si[iT+1], Se[iT+1] = Statistics(uf,A,𝗻ₖ,specie);
#         # end

#         # if savestats
#         #     flname = path*"/"*model_nm*"_statistics_t"*string(iT);
#         #     jldsave(flname, specie=specie, marg_labels=marg_labels, 
#         #     marg=marg[iT+1], 
#         #     E=𝔼[:,iT+1], Var=𝕍ar[:,iT+1],Sk=Sk[:,iT+1],
#         #      S=𝕊[iT+1], Si=Si[iT+1], Se=Se[iT+1], t=T[iT+1], T=T)
#         # end
#         # 𝕊[iT+1] = Entropy(uf);
#         # E = [sum(reshape(uf,𝗻ₖ...), dims=deleteat!(collect(1:length(𝗻ₖ)),i))[:] for i in 1:length(𝗻ₖ)];
#         # E = sum.(map(.*,E,collect.((:).(0,𝗻ₖ .- 1)))); 
#         # 𝔼[:,iT+1]=E;

#         # ProgressMeter.next!(pgres)
#     end

#     ub = copy(uf)[:]';

#     for iT in eachindex(T)[1:end-1]        
#         # prob = ODEProblem(g,ub, (T[end-iT+1],T[end-iT]));
#         prob = ODEProblem(g,ub, (T[iT],T[iT+1]));
#         sol = solve(prob, RK4();dt= .5/20, saveat= .5, adaptive=false);
#         ub = sol.u[end]/sum(sol.u[end]);

#         𝕊[end-iT] = Entropy(ub);
#         E = [sum(reshape(ub,𝗻ₖ...), dims=deleteat!(collect(1:length(𝗻ₖ)),i))[:] for i in 1:length(𝗻ₖ)];
#         E = sum.(map(.*,E,collect.((:).(0,𝗻ₖ .- 1)))); 
#         𝔼[:,end-iT]=E;

#         ProgressMeter.next!(pgres)
#     end

#     # 𝕊[end] = Entropy(uf);
#     # E = [sum(reshape(uf,𝗻ₖ...), dims=deleteat!(collect(1:length(𝗻ₖ)),i))[:] for i in 1:length(𝗻ₖ)];
#     # E = sum.(map(.*,E,collect.((:).(0,𝗻ₖ .- 1)))); 
#     # 𝔼[:,end]=E;


#     flname = path*"/"*model_nm*"_statistics";
#     jldsave(flname, E=𝔼, S=𝕊, T=T)
#     # return (marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊, Si, Se)
#     return 𝔼, 𝕊
# end