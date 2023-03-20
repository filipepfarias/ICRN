using NonlinearSolve, SteadyStateDiffEq, Sundials
using DifferentialEquations

model = "reactions/"*model_nm*".jl";
include(model);

𝐱₀ = [2e-7; 0; 5e-7; 0];

prob = ODEProblem(d𝐱dt!,𝐱₀,(T[1],T[end]));
sol = solve(prob,RK4(),saveat=.5, dt=.5, adaptive=false, progress=true);
nlprob = NonlinearProblem((𝐱,p,t) -> [d𝐱dt!(𝐱,p,t); sum(𝐱[2:end])-𝐱₀[3]; sum(𝐱[1:2])-𝐱₀[1]],sol(T[end]))
nlsol = solve(nlprob)
ssprob = SteadyStateProblem(prob);
sssol = solve(ssprob, dt = .01)

eₚ = [macroscopic_entropy_production(𝐱) for 𝐱 in sol.(T[2:end])];
𝑓d = [generalized_free_energy(𝐱,sol(T[end])) for 𝐱 in sol.(T[2:end])];
Eᵢₙ= [macroscopic_energy_input(𝐱,sol(T[end])) for 𝐱 in sol.(T[2:end])];

