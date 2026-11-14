function [V, f] = modes(K, M, n)
%MODES Lowest natural frequencies and mass-normalised mode shapes.
%   [V, f] = fmbd.modes(K, M, n) with f in Hz. K may be singular (free
%   systems); rigid-body modes then have f close to 0.

n = min(n, size(K, 1));
if size(K, 1) <= 200
    [V, L] = eig(full(K), full(M), 'vector');
else
    r = full(diag(K)./diag(M)); % shift below zero, scaled to the softest dofs
    r = r(r > 0 & isfinite(r)); % (keeps K - sigma*M well conditioned for sparse eigs)
    sigma = -1e-4;
    if ~isempty(r), sigma = -1e-4*min(r); end
    [V, L] = eigs(K, M, n, sigma);  L = diag(L);
end
[lam, i] = sort(real(L));
V = real(V(:, i(1:n)));  lam = lam(1:n);
V = V./sqrt(dot(V, M*V));
f = sqrt(max(lam, 0))/(2*pi);
end
