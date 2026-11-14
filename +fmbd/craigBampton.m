function [Kr, Mr, V, f] = craigBampton(K, M, b, n)
%CRAIGBAMPTON Craig-Bampton reduction with n fixed-interface modes.
%   Retains boundary dofs b; n = 0 gives rigid-body reduction. Returns reduced
%   matrices Kr, Mr, basis V, and fixed-interface frequencies f (Hz).
%   K(i,i), for non-boundary dofs i, must be nonsingular.
%
%   See: Craig, R. R., and Bampton, M. C. C., "Coupling of Substructures for
%   Dynamic Analyses," AIAA Journal, vol. 6, no. 7, 1968, pp. 1313-1319.

N = size(K, 1);
b = b(:);
i = setdiff((1:N)', b);
nb = numel(b);
Psi = -(decomposition(K(i,i), 'chol') \ full(K(i,b))); % constraint modes
if n > 0
    [Phi, f] = fmbd.modes(K(i,i), M(i,i), n); % fixed-interface modes
else
    Phi = zeros(numel(i), 0);  f = zeros(0, 1);
end
V = zeros(N, nb + n);
V(b, 1:nb) = eye(nb);
V(i, :) = [Psi, Phi];
Kr = V'*K*V;  Kr = (Kr + Kr')/2; % ensure that Kr is symmetric
Mr = V'*M*V;  Mr = (Mr + Mr')/2; % ensure that Mr is symmetric
end
