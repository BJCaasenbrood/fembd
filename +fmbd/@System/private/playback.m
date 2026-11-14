function playback(h, duration, frame)
% Call frame(t) with the elapsed wall-clock time until duration has passed.
t0 = tic;
while toc(t0) < duration && all(isvalid(h))
    frame(toc(t0));
    drawnow limitrate
end
if all(isvalid(h)) && isfinite(duration), frame(duration);  drawnow, end
end
