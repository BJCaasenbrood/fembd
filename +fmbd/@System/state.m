function x = state(s, q, nu, u)
% A one-step state struct as returned by simulate.
    if nargin < 4, u = zeros(numel(s.inSignals), 1); end
    x = struct('r', q.r, 'R', q.R, 'e', q.e, 'nu', nu, 'u', u);
end
