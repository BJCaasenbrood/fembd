function [f, V] = modes(s, x, n, opts)
%MODES Natural frequencies f (Hz) and mode shapes about a state at rest
%   (default: the initial configuration; e.g. from equilibrium):
%   [f, V] = s.modes(x, n). Unstable directions (e.g. an inverted
%   pendulum) give negative f. Option Input as in linearize. V (nDof x n) are velocity-like tangent
%   directions; s.plot(s.perturb(x, V(:,k)*a)) shows mode k.
    arguments
        s
        x struct = struct([])
        n (1,1) double = inf
        opts.Input double = []
    end
    lin = s.linearize(x, Input = opts.Input);
    [P, lam] = eig((lin.K + lin.K')/2, lin.M, 'vector');
    [lam, i] = sort(real(lam));  n = min(n, numel(lam));
    lam(abs(lam) < 1e-10*max(abs(lam))) = 0;
    f = sign(lam(1:n)).*sqrt(abs(lam(1:n)))/(2*pi);
    V = lin.T*real(P(:, i(1:n)));
end
