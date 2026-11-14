function h = plotMode(s, x, v, opts)
%PLOTMODE Draw a mode shape: s.plotMode(x0, V(:,k)) with x0 and V from
%   modes, displaced so that the largest node displacement is Amplitude
%   (default 0.05 m). Option Title. Returns one handle per body.
    arguments
        s
        x struct
        v (:,1) double
        opts.Amplitude (1,1) double {mustBePositive} = 0.05
        opts.Title = []
    end
    [x, a] = s.modeScale(x, v, opts.Amplitude);
    h = s.plot(s.perturb(x, a*v));
    if ~isempty(opts.Title), title(opts.Title); end
    if nargout == 0, clear h, end
end
