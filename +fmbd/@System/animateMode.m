function animateMode(s, x, v, opts)
%ANIMATEMODE Animate a mode shape: s.animateMode(x0, V(:,k)) with x0
%   and V from modes, oscillating about x0 in slow motion (vibration
%   modes are too fast for real time). The phase follows the wall clock.
%   Options: Amplitude (largest node displacement, default 0.05 m),
%   Frequency (Hz, shown with the slow-motion factor), Rate (cycles per
%   second, default 1), Duration (s, default 3), Title.
    arguments
        s
        x struct
        v (:,1) double
        opts.Amplitude (1,1) double {mustBePositive} = 0.05
        opts.Frequency double = []
        opts.Rate (1,1) double {mustBePositive} = 1
        opts.Duration (1,1) double {mustBePositive} = 3
        opts.Title = []
    end
    [x, a, V0] = s.modeScale(x, v, opts.Amplitude);  nb = numel(s.Bodies);
    cla;
    h = drawBodies(s.Bodies, V0, num2cell(lines(nb), 2)');
    axis manual
    ttl = char(opts.Title);
    if ~isempty(opts.Frequency)
        ft = sprintf('%.1f Hz (%.0fx slow motion)', opts.Frequency, opts.Frequency/opts.Rate);
        if opts.Frequency <= 0, ft = 'rigid-body / mechanism mode'; end
        if isempty(ttl), ttl = ft; else, ttl = [ttl ': ' ft]; end
    end
    if ~isempty(ttl), title(ttl); end
    playback(h, opts.Duration, @(c) setVertices(h, arrayfun(@(i) s.vertices(i, ...
        s.perturb(x, a*sin(2*pi*opts.Rate*c)*v), 1, 1), 1:nb, 'UniformOutput', false)));
end
