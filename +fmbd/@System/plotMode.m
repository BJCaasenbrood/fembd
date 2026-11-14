function h = plotMode(s, x, v, opts)
%PLOTMODE Draw a mode shape.
%   s.plotMode(x0, v) with x0 and v = V(:,k) from modes.
%   Options: Amplitude (largest node displacement, default 0.05), Title.
    arguments
        s
        x struct
        v (:,1) double
        opts.Amplitude (1,1) double {mustBePositive} = 0.05
        opts.Title = []
    end
    [x, a] = modeScale(s, x, v, opts.Amplitude);
    h = s.plot(s.perturb(x, a*v));
    if ~isempty(opts.Title), title(opts.Title); end
    if nargout == 0, clear h, end
end
