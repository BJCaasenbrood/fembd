function [v, Jv] = coord(s, kind, idx, q)
% Coordinates and Jacobian of a joint (its free coordinates), a
% reference point ([displacement; rotation vector]) or a node set
% (average displacement), all zero in the defining configuration.
    switch kind
        case 'joint'
            [c, Jc] = s.jointKin(idx, q);  f = s.Jn(idx).free;
            v = c(f);  Jv = Jc(f, :);
        case 'ref'
            R = s.Rf(idx);
            [ra, Ra, Ja, Jw] = s.anchorKin(R.anchor, q);
            phi = logso3(Ra);
            v = [ra + Ra*R.s - R.x0; phi];
            Jv = [Ja - Ra*hat(R.s)*Jw; dexpInv(phi)*Jw];
        otherwise                            % node set
            Ns = s.Nd(idx);  B = s.Bd(Ns.body);  c = B.cols;
            Rb = q.R(:,:,Ns.body);  ub = Ns.u0 + Ns.Phi*q.e(B.ecols);
            v = q.r(:, Ns.body) + Rb*ub - Ns.x0;
            Jv = zeros(3, s.N);
            Jv(:, c(1:3)) = eye(3);  Jv(:, c(4:6)) = -Rb*hat(ub);  Jv(:, c(7:end)) = Rb*Ns.Phi;
    end
end
