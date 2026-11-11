function playback(h, duration, frame)
%PLAYBACK Call frame(s) with the elapsed wall-clock time s (seconds) until
%   duration has passed or a graphics object in h is closed, then once more
%   at the end, so the speed does not depend on the rendering speed.
t0 = tic;
while toc(t0) < duration && all(isvalid(h))
    frame(toc(t0));
    drawnow limitrate
end
if all(isvalid(h)) && isfinite(duration), frame(duration);  drawnow, end
end
