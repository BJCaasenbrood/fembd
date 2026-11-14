function [V, f] = modes(K, M, n)
%MODES Lowest n natural frequencies f (Hz) and mass-normalised mode shapes V.
%   Works for free (singular K) systems: rigid-body and mechanism modes
%   appear as f ~ 0. Small systems are solved densely.

n = min(n, size(K, 1));
if size(K, 1) <= 200
    [V, L] = eig(full(K), full(M), 'vector');
else
    r = full(diag(K)./diag(M)); % shift below zero, scaled to the softest
    r = r(r > 0 & isfinite(r)); % dofs (keeps K - sigma*M well conditioned for sparse eigs)
    sigma = -1e-4;
    if ~isempty(r), sigma = -1e-4*min(r); end
    [V, L] = eigs(K, M, n, sigma);  L = diag(L);
end
[lam, i] = sort(real(L));
V = real(V(:, i(1:n)));  lam = lam(1:n);
V = V./sqrt(dot(V, M*V));
f = sqrt(max(lam, 0))/(2*pi);
end
