function add_probability_condition(A)
    return [A[1:end-1,:]; ones(1,size(A,2))]
end

function get_prob_eq(A,p)
    y = 0.0*p;
    y[end] = 1.0;
    return add_probability_condition(A)\y
end

function Gillespie(x0,K,rs,tspan)
    𝛎 = netstoichmat(rs);
    Re = substoichmat(rs);
    T = tspan[end];
    t = Vector{Float32}();
    tt = tspan[1];
    xx = x0;
    x = Vector{Any}();
    while tt <= T
        𝛂ₘ = prod(binomial.(xx,Re),dims=1);
        𝛂 = (tt .|> K) .* 𝛂ₘ;

        α₀ = sum(𝛂);
        α₀ == 0.0 ? break : nothing
        tt += log(1 / rand()) / α₀;
        xx += 𝛎[:,sample(Weights(vec(𝛂/α₀)))];
        push!(t,tt);
        push!(x,xx);
    end
    return (t,x)
end