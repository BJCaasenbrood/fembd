function [v, Jv] = coord(s, kind, idx, q)
% Coordinates and Jacobian of a joint, reference point or node set.
    switch kind
        case 'joint'
            [c, Jc] = jointKin(s, idx, q);  f = s.jointFrames(idx).free;
            v = c(f);  Jv = Jc(f, :);
        case 'ref'
            R = s.refFrames(idx);
            [ra, Ra, Ja, Jw] = anchorKin(s, R.anchor, q);
            phi = logso3(Ra);
            v = [ra + Ra*R.s - R.x0; phi];
            Jv = [Ja - Ra*hat(R.s)*Jw; dexpInv(phi)*Jw];
        otherwise % node set
            Ns = s.nodeSets(idx);
            B = s.units(Ns.body);
            c = B.cols;
            Rb = q.R(:,:,Ns.body);  ub = Ns.u0 + Ns.Phi*q.e(B.ecols);
            v = q.r(:, Ns.body) + Rb*ub - Ns.x0;
            Jv = zeros(3, s.nv);
            Jv(:, c(1:3)) = eye(3);  Jv(:, c(4:6)) = -Rb*hat(ub);  Jv(:, c(7:end)) = Rb*Ns.Phi;
    end
end
