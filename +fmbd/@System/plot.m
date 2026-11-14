function h = plot(s, x, k, opts)
%PLOT Draw the system.
%   h = s.plot() draws the initial configuration, s.plot(x, k) step k of a
%   simulation (default: the last). Option Scale magnifies the deformation.
    arguments
        s
        x = []
        k = []
        opts.Scale (1,1) double = 1
    end
    s.build();
    if isempty(x), x = s.state(s.q0, s.nu0); end
    if isempty(k), k = size(x.nu, 2); end
    nb = numel(s.Bodies);
    V = arrayfun(@(i) s.vertices(i, x, k, opts.Scale), 1:nb, 'UniformOutput', false);
    h = drawBodies(s.Bodies, V, num2cell(lines(nb), 2)');
    if nargout == 0, clear h, end
end
