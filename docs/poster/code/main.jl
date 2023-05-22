using LinearAlgebra, DifferentialEquations
using CairoMakie
include("utils.jl")

# k₋ = 1/31e-3;
k₋ = .1 + 5e-3;
# γ = .16;
γ = 0.0;
tₘ = 40e-0;
λ(t) = γ*sin(2π*t/tₘ);
# β = 1/26e-3;
β = .5;

A(t) = [ -β*(1+λ(t)) k₋;
          β*(1+λ(t)) -k₋];

function update_func(B_,u,p,t)
    B_[:] .= A(t)[:];
    # B_[1,1] = A(t)[1,1];
    # B_[2,1] = A(t)[2,1];
    # B_[1,2] = A(t)[1,2];
    # B_[2,2] = A(t)[2,2];
end

B_ = DiffEqArrayOperator(ones(2,2), update_func = update_func)

# prob = ODEProblem(B_, ones(2)/2, (0.0, 800e-3))
prob = ODEProblem(B_, ones(2)/2, (0.0, 25))
# sol = solve(prob,MagnusGL6(),dt = 1e-3)
sol = solve(prob,MagnusGL6(),dt = 5e-1)


# Entropy balace
Ṡ_sys(p,Q) =  .5*sum(x -> (isnan(x) && x<=0.0) ? 0 : x ,(diagm(p)Q'-Q*diagm(p)) .* log.(p ./ p') );
Ṡ_bath(p,Q) = .5*sum(x -> (isnan(x) && x<=0.0) ? 0 : x ,(diagm(p)Q'-Q*diagm(p)) .* (log.(Q' ./ Q)) );
Ṡ_tot(p,Q) = .5*sum(x -> (isnan(x) && x<=0.0) ? 0 : x ,(diagm(p)Q'-Q*diagm(p)) .* (log.((diagm(p)Q') ./ (Q*diagm(p)))) );
Ḟ(p,Q) = sum(x -> (isnan(x) && x<=0.0) ? 0 : x ,Q'p .* log.( p ./ get_prob_eq(Q,p) ) );
F(p,Q) = sum(x -> (isnan(x) && x<=0.0) ? 0 : x ,p .* (log.( p ./ get_prob_eq(Q,p) )));
# ẇ(p,Q) = -.5*sum(x -> (isnan(x) && x<=0.0) ? 0 : x ,(diagm(p)Q'-Q*diagm(p)) .* (log.((diagm(get_prob_eq(Q,p))Q') ./ (Q*diagm(get_prob_eq(Q,p))))) );

vec_Ṡ_sys = [Ṡ_sys(sol(t),A(t)) for t in sol.t];
vec_Ṡ_bath = [Ṡ_bath(sol(t),A(t)) for t in sol.t];
vec_Ṡ_tot = [Ṡ_tot(sol(t),A(t)) for t in sol.t];
vec_Ḟ = [Ḟ(sol(t),A(t)) for t in sol.t];
vec_F = [F(sol(t),A(t)) for t in sol.t];
# vec_ẇ = [ẇ(sol(t),A(t)) for t in sol.t];


f = Figure()
ax1 = Axis(f[1, 1], title = "Entropy rates")

# labels=["Ṡ_sys" "Ṡ_bath" "Ṡ_tot" "F" "Ḟ"]
labels=["Ṡ_sys" "Ṡ_bath" "Ṡ_tot" "Ḟ"]
# for (i,y) in enumerate([vec_Ṡ_sys,vec_Ṡ_bath,vec_Ṡ_tot,vec_F,vec_Ḟ])
for (i,y) in enumerate([vec_Ṡ_sys,vec_Ṡ_bath,vec_Ṡ_tot,vec_Ḟ])
    # if i == 5
    if false
        scatterlines!(ax1,sol.t,y, label = labels[i],markersize = 10, marker = :cross, strokecolor = :white, strokewidth = .5)
    else
        scatterlines!(ax1,sol.t,y, label = labels[i], strokecolor = :white, strokewidth = .5)
    end
end
axislegend(ax1)

labels=["p₀" "p₁" "pᵉ₀" "pᵉ₁"]
ax2 = Axis(f[2, 1], title="Time-evolution of the Probabilities")
for (i,y) in enumerate((hcat(sol.u...)'[:,1],hcat(sol.u...)'[:,2],hcat([get_prob_eq(A(t),sol(t)) for t in sol.t]...)'[:,1], hcat([get_prob_eq(A(t),sol(t)) for t in sol.t]...)'[:,2]))
    if i <= 2
        scatterlines!(ax2,sol.t,y, label = labels[i], strokecolor = :white, strokewidth = .5)
    else
        lines!(ax2,sol.t,y, label = labels[i],linestyle = :dash)
    end
end
axislegend(ax2)
f