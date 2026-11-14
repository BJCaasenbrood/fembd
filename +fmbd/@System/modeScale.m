function [x, a, V0] = modeScale(s, x, v, amplitude)
% State x (default: initial), the factor a that makes the largest node
% displacement of a*v equal to amplitude, and the node positions V0.
    s.build();  nb = numel(s.Bodies);
    if isempty(x), x = s.state(s.q0, s.nu0); end
    V0 = arrayfun(@(i) s.vertices(i, x, size(x.nu, 2), 1), 1:nb, 'UniformOutput', false);
    e = 1e-6;  x1 = s.perturb(x, e*v);         % node motion per unit amplitude
    D = arrayfun(@(i) (s.vertices(i, x1, 1, 1) - V0{i})/e, 1:nb, 'UniformOutput', false);
    a = amplitude/max([cellfun(@(d) max(vecnorm(d, 2, 2)), D), eps]);
end
