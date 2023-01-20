using CairoMakie

fig = Figure(resolution = (300,300));

fig, ax2, sp = series(sol_Det[3],sol_CME[3], labels=sol_CME[1][[1, 5, 8, 10]] .* "-CME",linestyle = :dot);
series!(ax2,sol_Det[3],sol_SSA[3], labels=sol_SSA[1][[1, 5, 8, 10]] .* "-SSA");
series!(ax2,sol_Det[3],sol_Det[2], labels=sol_Det[1] .* "-Det",markersize=6);
axislegend(ax2);
fig

mkdir(path*"/plots")

save(path*"/plots/"*model_nm*"_mean_evol.pdf", fig, pt_per_unit = 2)

fig = Figure(resolution = (300,300));

fig, ax2, sp = series(sol_Det[3],sol_CME[6], labels=["Entropy_CME"],linestyle = :dot);
series!(ax2,sol_Det[3],sol_SSA[6], labels=["Entropy_SSA"]);
axislegend(ax2);
fig

save(path*"/plots/"*model_nm*"_mean_evol.pdf", fig, pt_per_unit = 2)
