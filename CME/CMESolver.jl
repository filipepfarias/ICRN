using DifferentialEquations

include("src/CME.jl");
include("reactions/MichaelisMenten.jl");

T = 0.0:.5:50.0;
f(u,p,t) = A*u ;

begin
    u0 = p₀[:];
    p = Vector{Vector{Float64}}([u0]);
    for iT in eachindex(T)[1:end-1]
        global u0
        prob = ODEProblem(f,u0, (T[iT],T[iT+1]));
        sol = solve(prob, RK4(),saveat=T[iT+1]);
        sol.u[end][sol.u[end] .< 0] .= 0;
        append!(p,[sol.u[end]]);
        u0 = sol.u[end]/sum(sol.u[end]);
    end
end



