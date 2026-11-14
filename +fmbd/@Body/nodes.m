function idx = nodes(b, sel)
%NODES Indices of nodes for which sel(X) is true (X is n x 3), or the
%   node indices sel themselves if sel is numeric.
    if isnumeric(sel), idx = sel(:); else, idx = find(sel(b.X)); end
end
