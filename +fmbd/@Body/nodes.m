function idx = nodes(b, sel)
%NODES Node indices selected by sel(X), or sel itself if numeric.
    if isnumeric(sel), idx = sel(:); else, idx = find(sel(b.X)); end
end
