function q = retract(s, q, dq)
% q o exp(dq): positions add (world), rotations update on the right
% (body), modal coordinates add.
    for i = 1:numel(s.units)
        c = s.units(i).cols;
        q.r(:, i) = q.r(:, i) + dq(c(1:3));
        q.R(:,:,i) = q.R(:,:,i)*expso3(dq(c(4:6)));
        q.e(s.units(i).ecols) = q.e(s.units(i).ecols) + dq(c(7:end));
    end
end
