function G = gapOperator(t, nn, masterOnly)
% Sparse gap operator G*u = u_slave - sum_j N_j u_master,j of a tie (rows
% at the slave dofs); with masterOnly, the master part only.
if nargin < 3, masterOnly = false; end
d = reshape(1:3, 1, 1, 3);  ns = numel(t.snodes);
if masterOnly
    rows = 3*(t.snodes(:) - 1) + d + 0*t.mnodes;  cols = 3*t.mnodes - 3 + d;
    vals = -t.N + 0*d;
    G = sparse(rows(:), cols(:), vals(:), 3*nn, 3*nn);
else
    nodes = [t.snodes(:), t.mnodes];
    rows = 3*(0:ns - 1)' + d + 0*nodes;  cols = 3*nodes - 3 + d;
    vals = [ones(ns, 1), -t.N] + 0*d;
    G = sparse(rows(:), cols(:), vals(:), 3*ns, 3*nn);
end
end
