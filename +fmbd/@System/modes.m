function [f, V] = modes(s, x, n, opts)
%MODES Natural frequencies and mode shapes about a state at rest.
%   [f, V] = s.modes(x, n) returns n frequencies (Hz) and mode shapes.
%   Unstable directions have negative f. Option Input as in linearize.
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
