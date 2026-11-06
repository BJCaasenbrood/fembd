function [Kr, Mr, V, f] = craigBampton(K, M, b, n)
%CRAIGBAMPTON Craig-Bampton reduction with boundary dofs b and n modes.
%   [Kr, Mr, V, f] = fem.craigBampton(K, M, b, n) keeps the dofs b (indices
%   into K, e.g. from Assembly.qdof) as physical coordinates and adds n
%   fixed-interface modes. The reduced coordinates are x = [u_b; eta] with
%   u ~ V*x, Kr = V'*K*V and Mr = V'*M*V. f are the fixed-interface
%   frequencies (Hz). n = 0 gives Guyan (static) reduction.
%
%   K(i,i), with i the non-boundary dofs, must be non-singular: b (together
%   with any constraints already eliminated from K) must stop rigid motion.

N = size(K, 1);  b = b(:);  i = setdiff((1:N)', b);  nb = numel(b);
Psi = -(decomposition(K(i,i), 'chol') \ full(K(i,b)));    % constraint modes
if n > 0
    [Phi, f] = fem.modes(K(i,i), M(i,i), n);              % fixed-interface modes
else
    Phi = zeros(numel(i), 0);  f = zeros(0, 1);
end
V = zeros(N, nb + n);
V(b, 1:nb) = eye(nb);
V(i, :) = [Psi, Phi];
Kr = V'*K*V;  Kr = (Kr + Kr')/2;
Mr = V'*M*V;  Mr = (Mr + Mr')/2;
end
