function lotus(g,p::SparseVector{<:Number,Int})
    return sum(g.(p.nzind) .* p.nzval);
end

function sparselog(p::AbstractSparseArray)
    nzp = nonzeros(p);
    nzp .= log.(nzp);
return p
end

function marginal(p::SparseVector{<:Number,Int},𝗻ₖ,dim)
    # Tensor contraction fashion
    𝓘 = getindex.(getfield.(CartesianIndices(𝗻ₖ)[p.nzind],:I),dim);
    pmarg = spzeros(𝗻ₖ[dim]);
    for i in 1:𝗻ₖ[dim];
        pmarg[i] = sum(p.nzval[𝓘 .== i])
    end
    return pmarg
end

function mean(p::SparseVector{<:Number,Int},𝗻ₖ,dim)
    pmarg = marginal(p,𝗻ₖ,dim)
    return lotus(x -> identity(x - 1),pmarg)
end

# function variance(p::SparseVector{<:Number,Int},𝗻ₖ,dim)
#     pmarg = marginal(p,𝗻ₖ,dim)
#     return lotus(x -> identity(x - 1),pmarg)
# end

function entropy(p::SparseVector{<:Number,Int})
    nzlog(x) = x <= 0 ? 0 : log(x);
    g(i) = nzlog(p[i]);
    return -lotus(g,p)
end

function d_entropy(p::SparseVector{<:Number,Int},Q::SparseMatrixCSC{<:Real,<:Int64})
    return -sum(sparselog(p) .* (Q*p))
end

function KLdivergence(p::SparseVector{<:Number,Int},q::SparseVector{<:Number,Int})
    nzlog(x) = x <= 0 ? 0 : log(x);
    g(x) = nzlog(p[x]/q[x]);
    return lotus(g,p)
end

function entropy_flow(p::SparseVector{<:Number,Int},Q::SparseMatrixCSC{<:Real,<:Int64})
    Q = Q - spdiagm(diag(Q));
    j = spdiagm(p)Q';
    J = j - j';
    logQ = copy(Q');
    nzlogQ = nonzeros(logQ);
    nzlogQ .= log.(nzlogQ);
    return .5*sum(J .* (logQ - logQ'));
end

function entropy_production(p::SparseVector{<:Number,Int},Q::SparseMatrixCSC{<:Real,<:Int64})
    Q = Q - spdiagm(diag(Q));
    j = spdiagm(p)Q';
    J = j - j';
    logj = copy(j);
    nzlogj = nonzeros(logj);
    nzlogj .= log.(nzlogj);
    return .5*sum(J .* (logj - logj'));
end

function free_energy(p::SparseVector{<:Number,Int},pss::SparseVector{<:Number,Int})
    return KLdivergence(p,pss)
end

function energy_input_rate(pss::SparseVector{<:Number,Int},p::SparseVector{<:Number,Int},Q::SparseMatrixCSC{<:Real,<:Int64})
    Q = Q - spdiagm(diag(Q));
    j = spdiagm(p)Q';
    J = j - j';
    j = spdiagm(pss)Q';
    logj = copy(j);
    nzlogj = nonzeros(logj);
    nzlogj .= log.(nzlogj);
    return .5*sum(J .* (logj - logj'));
end

function time_reversed_entropy_production(p::SparseVector{<:Number,Int},Q::SparseMatrixCSC{<:Real,<:Int64})
    Q = Q - spdiagm(diag(Q));
    J = (Q * spdiagm(p));
    logJ = copy(J);
    nzJ = nonzeros(logJ);
    nzJ .= log.(nzJ);
    return sum((J - J') .* (logJ - logJ'))
end

function time_reversed_entropy_flow(p::SparseVector{<:Number,Int},Q::SparseMatrixCSC{<:Real,<:Int64})
    Q = Q - spdiagm(diag(Q));
    logQ = copy(Q);
    nzQ = nonzeros(logQ);
    nzQ .= log.(nzQ);
    J = (Q * spdiagm(p));
    return sum((J - J') .* (logQ - logQ'))
end


# function Entropy(p)
#     ip = p .!= 0.0
#     𝕊 = -p[ip] .* log.(p[ip])
#     𝕊 = sum(𝕊);

#     return 𝕊
# end

# function dEntropy(p,A)
#     Q = A - spdiagm(diag(A));
#     logA = copy(Q);
#     nzlogA = nonzeros(logA); nzlogA .= log.(nonzeros(logA));

#     logAp = Q*dropzeros(spdiagm(p));
#     nzlogAp = nonzeros(logAp); nzlogAp .= log.(nonzeros(logAp));

#     J = Q*dropzeros(spdiagm(p)) - (Q*dropzeros(spdiagm(p)))';
#     X = logAp - logAp';

#     Se = .5 * sum(J .* (logA - logA')); # h_ex
#     Si = .5 * sum( J .* X ); # e_p

#     return Si, Se
# end

# # function CMEMutualInformation(Xₖ₋₁,Xₖ,A,dt)
    
# #     ℐ = 0.0;

# #     for i in eachindex(Xₖ), j in eachindex(Xₖ₋₁)
# #         expA = (I[i,j]+A[i,j]*dt);
# #         logXₖₖ₋₁ = log(expA*Xₖ₋₁[j]);
# #         logXₖ₋₁  = log(Xₖ₋₁[j]);
# #         logXₖ   = log(Xₖ[i]);
        
# #         if !any(isinf.([logXₖₖ₋₁, logXₖ₋₁, logXₖ]))
# #             ℐ += expA*Xₖ₋₁[j] * (logXₖₖ₋₁ - logXₖ₋₁ - logXₖ)
# #         end
# #     end

# #     return ℐ
# # end

# function Mean(𝗻ₖ,marg)
#     𝔼 = zeros(length(𝗻ₖ),1);

#     for i in eachindex(𝗻ₖ)
#         ℕ = 0:(𝗻ₖ[i]-1);
#         𝔼[i,1] = sum( ℕ .* marg[i] );
#     end
#     return 𝔼
# end

# function Variance(𝗻ₖ,𝔼,marg)
#     𝕍ar = zeros(length(𝗻ₖ),1);

#     for i in eachindex(𝗻ₖ)
#         ℕ = (1:𝗻ₖ[i]) .- 1;
#         𝕍ar[i,1] = sum( (ℕ.-𝔼[i,1]).^2 .* marg[i] );
#     end
#     return 𝕍ar
    
# end

# function Skewness(𝗻ₖ,𝔼,𝕍ar,marg)
#     Sk = zeros(length(𝗻ₖ),1);

#     for i in eachindex(𝗻ₖ)
#         ℕ = (1:𝗻ₖ[i]) .- 1;
#         Sk[i,1] = 𝕍ar[i,1] == 0.0 ? 0.0 : sum(((ℕ.-𝔼[i,1])./√𝕍ar[i,1]).^3 .* marg[i] );
#     end
#     return Sk
# end

# function Marginals(𝗻ₖ,𝓅,specie)
#     idn = Int(length(𝗻ₖ)*(length(𝗻ₖ)+1)/2);
#     marg = Vector{Any}(undef,idn);
#     idm = 1;
#     marg_labels = [];
#     # 𝓅ₙ = sum(𝓅);
#     𝓅ₙ = 1;

#     for i in eachindex(𝗻ₖ), j in eachindex(𝗻ₖ)
#         if j > i
#             # ℕxℕ = ((1:𝗻ₖ[i]) .- 1)*((1:𝗻ₖ[i]) .- 1)';
#             d = deleteat!(collect(eachindex(𝗻ₖ)), [i j])
#             marg[idm] = reshape(sum(𝓅,dims=d)./𝓅ₙ ,𝗻ₖ[i],𝗻ₖ[j])

#             # ℝ[i,1,j] = sum( ℕxℕ .* marg[idm] );

#             flsuffix = specie[i]*"_x_"*specie[j];
#             append!(marg_labels,[flsuffix]);
#             idm += 1;
#         elseif i == j
#             ind = collect(1:length(𝗻ₖ))
#             marg[idm] = sum(𝓅,dims=deleteat!(ind,i))[:] ./𝓅ₙ ;

#             flsuffix = specie[i];
#             append!(marg_labels,[flsuffix]);
#             idm += 1;
#         end
#     end

#     return marg, marg_labels
# end

# function Statistics(p,A,𝗻ₖ,specie)
#     𝓅  = reshape(p,𝗻ₖ...,);

#     marg, marg_labels = Marginals(𝗻ₖ,𝓅,specie);
#     i_marg = [1; 1 .+ cumsum(length(𝗻ₖ):-1:2)]
#     𝔼 = Mean(𝗻ₖ,marg[i_marg])
#     𝕍ar = Variance(𝗻ₖ,𝔼,marg[i_marg]);
#     Sk = Skewness(𝗻ₖ,𝔼,𝕍ar,marg[i_marg]);

#     𝕊 = Entropy(p);

#     Si, Se = dEntropy(p,A);
#     return marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊, Si, Se
# end

# function Statistics(p,𝗻ₖ,specie)
#     𝓅  = reshape(p,𝗻ₖ...,);

#     marg, marg_labels = Marginals(𝗻ₖ,𝓅,specie);
#     i_marg = [1; 1 .+ cumsum(length(𝗻ₖ):-1:2)]
#     𝔼 = Mean(𝗻ₖ,marg[i_marg])
#     𝕍ar = Variance(𝗻ₖ,𝔼,marg[i_marg]);
#     Sk = Skewness(𝗻ₖ,𝔼,𝕍ar,marg[i_marg]);

#     𝕊 = Entropy(p);

#     return marg_labels, marg, 𝔼, 𝕍ar, Sk, 𝕊
# end

# using FileIO, JLD2

# function saveStatistics(path, model_nm)

#     include(path*"/model.jl");
#     # A = CMEOperator(𝛎,Re,K,𝗻ₖ);
    
#     marg_labels = [];
#     marg = Vector{Any}(undef,length(T));
#     𝔼 = zeros(length(𝗻ₖ),length(T));
#     𝕍ar = zeros(length(𝗻ₖ),length(T));
#     Sk = zeros(length(𝗻ₖ),length(T));
#     𝕊 = zeros(1,length(T));
#     # Si = zeros(1,length(T));
#     # Se = zeros(1,length(T));

#     for iT in eachindex(T)        
#         data = jldopen(path*"/"*model_nm*"_t"*string(iT-1));
#         p = data["p"];

#         # marg_labels, marg[iT], 𝔼[:,iT], 𝕍ar[:,iT], Sk[:,iT], 𝕊[iT], Si[iT], Se[iT] = Statistics(p,A,𝗻ₖ,specie);
#         marg[i], _ = Marginals(𝗻ₖ,P,specie);
#         i_marg = [1; 1 .+ cumsum(length(𝗻ₖ):-1:2)]
#         𝔼[:,i] = Mean(𝗻ₖ,marg[i][i_marg])
#         𝕍ar[:,i] = Variance(𝗻ₖ,𝔼,marg[i][i_marg]);
#         Sk[:,i] = Skewness(𝗻ₖ,𝔼,𝕍ar,marg[i][i_marg]);

#         𝕊[i] = Entropy(P);
#     end

#     flname = path*"/"*model_nm*"_statistics";
#     # jldsave(flname, specie=specie, marg_labels=marg_labels, marg=marg, E=𝔼, Var=𝕍ar,Sk=Sk,
#     #          S=𝕊, Si=Si, Se=Se, t=T, T=T)
#     jldsave(flname, specie=specie, marg_labels=marg_labels, marg=marg, E=𝔼, Var=𝕍ar,Sk=Sk,
#             S=𝕊, t=T, T=T)
# end