function animate(s, t, x, opts)
%ANIMATE Animate a simulation: s.animate(t, x) with t, x from simulate.
%   Options: Rate (default 1) = simulated seconds per wall-clock
%   second, so 1 is real time and 0.1 slow motion; Scale (default 1)
%   magnifies the elastic deformation; Title.
    arguments
        s
        t (:,1) double
        x struct
        opts.Rate (1,1) double {mustBePositive} = 1
        opts.Scale (1,1) double = 1
        opts.Title (1,:) char = ''
    end
    s.build();
    nb = numel(s.Bodies);
    V = @(k) arrayfun(@(i) s.vertices(i, x, k, opts.Scale), 1:nb, 'UniformOutput', false);
    cla;  h = s.plot(x, 1, Scale = opts.Scale);
    X = reshape(permute(x.r, [1 3 2]), 3, [])';           % all frame positions
    span = max([arrayfun(@(B) max([vecnorm(B.Xl, 2, 2); 0]), s.Bd), 1e-3]);
    axis(reshape([min(X, [], 1) - span; max(X, [], 1) + span], 1, []));
    playback(h, (t(end) - t(1))/opts.Rate, @(c) frame(c));
    function frame(c)
        k = find(t <= t(1) + opts.Rate*c, 1, 'last');
        setVertices(h, V(k));
        title(sprintf('%s  t = %.3f s', opts.Title, t(k)));
    end
end
