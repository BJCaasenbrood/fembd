function p = logso3(R)
%LOGSO3 Rotation vector of a rotation matrix (angle in [0, pi]).
s = [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)]/2;
c = (trace(R) - 1)/2;
sn = norm(s);
th = atan2(sn, c);
if th < 1e-8
    p = s;
elseif pi - th > 1e-6
    p = th/sn*s;
else % near pi: axis from R + I
    B = (R + eye(3))/2;  [~, i] = max(diag(B));
    n = B(:, i)/sqrt(B(i, i));
    if n'*s < 0, n = -n; end
    p = th*n;
end
end
