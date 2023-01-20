function concentration2mol(K,Re,V,binv=1)
    nₐ  = 6.022e23; # Avogadro number
    c = factorial.(Re);
    # c[Re .== 0] = 1;
    c = prod(c, dims = 2);

    i0 = map(x -> all(x .== 0),eachrow(Re)); # 0th order
    c[i] .= K[i] * (nₐ * V)^(binv);

    i1 = map(x -> sum(x) == 1,eachrow(Re)); # 1st order
    c[i] .= K[i];

    i0i1 = i0 .* i1; # Higher order
    c[i0i1] .= K[i0i1] * (r[i0i1]).^(binv) / (nₐ * V)^(binv)

    return c
end

function mol2concentration(K,Re,V,binv=-1)
    K = concentration2mol(K,Re,V,binv)
end