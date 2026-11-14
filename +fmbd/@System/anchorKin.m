function [ra, Ra, Jv, Jw] = anchorKin(s, k, q)
% Position ra and rotation Ra of anchor k (0 = ground) with the
% Jacobians of its velocity (world) and angular velocity (anchor frame).
    Jv = zeros(3, s.N);  Jw = zeros(3, s.N);
    if k == 0, ra = zeros(3, 1);  Ra = eye(3);  return, end
    A = s.An(k);  B = s.Bd(A.body);  c = B.cols;
    r = q.r(:, A.body);  R = q.R(:,:,A.body);
    Jv(:, c(1:3)) = eye(3);
    if A.frame
        ra = r;  Ra = R;  Jw(:, c(4:6)) = eye(3);
        return
    end
    eta = q.e(B.ecols);  ut = A.x0 + A.Psi(1:3,:)*eta;  psi = A.Psi(4:6,:)*eta;  E = expso3(psi);
    ra = r + R*ut;  Ra = R*E;
    Jv(:, c(4:6)) = -R*hat(ut);  Jv(:, c(7:end)) = R*A.Psi(1:3,:);
    Jw(:, c(4:6)) = E';          Jw(:, c(7:end)) = dexp(psi)*A.Psi(4:6,:);
end
