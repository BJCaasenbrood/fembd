function [F, Kt, Ct, Mq] = forces(s, q, nu, u, tangents)
% Generalised forces F (gravity, inertia, elastic, springs, dampers,
% inputs), approximate tangents Kt = -dF/dq, Ct = -dF/dnu (only when
% tangents is true, the default), and M(q).
    if nargin < 5, tangents = true; end
    n = s.N;  F = zeros(n, 1);  Kt = zeros(n);  Ct = zeros(n);  Mq = zeros(n);
    g = s.Gravity(:);
    for i = 1:numel(s.Bd)
        B = s.Bd(i);  c = B.cols;  w = nu(c(4:6));  R = q.R(:,:,i);
        if B.elastic
            eta = q.e(B.ecols);
            [Mb, Fb] = elasticTerms(B, R, eta, nu(c), g);
            Mq(c, c) = Mb;  F(c) = Fb;
            if tangents                      % velocity tangent: exact (quadratic in nu)
                Cb = zeros(numel(c), numel(c));
                for k = 1:numel(c)
                    dv = zeros(numel(c), 1);  dv(k) = 1;
                    [~, Fp] = elasticTerms(B, R, eta, nu(c) + dv, g);
                    [~, Fm] = elasticTerms(B, R, eta, nu(c) - dv, g);
                    Cb(:, k) = (Fm - Fp)/2;
                end
                Ct(c, c) = Cb;
                Kt(c(7:end), c(7:end)) = B.Kc;
            end
        else
            Mq(c, c) = blkdiag(B.m*eye(3), B.J);
            F(c(1:3)) = B.m*g;
            F(c(4:6)) = -cross(w, B.J*w);
            Ct(c(4:6), c(4:6)) = hat(w)*B.J - hat(B.J*w);
        end
    end
    for j = 1:numel(s.Jn)
        J = s.Jn(j);  f = J.free & (J.k' ~= 0 | J.c' ~= 0);
        if ~any(f), continue, end
        [c, Jc] = s.jointKin(j, q);  c = c(f);  Jc = Jc(f, :);  k = J.k(f);  d = J.c(f);
        F = F - Jc'*(k.*c + d.*(Jc*nu));
        Kt = Kt + Jc'*(k.*Jc);  Ct = Ct + Jc'*(d.*Jc);
    end
    for T = s.Tr(arrayfun(@(T) ~isempty(T.k), s.Tr))
        [d, Jd] = s.stretch(T, q);
        F = F - Jd'*(T.k*d + T.c*(Jd*nu));
        Kt = Kt + T.k*(Jd'*Jd);  Ct = Ct + T.c*(Jd'*Jd);
    end
    for i = 1:numel(s.In)
        sg = s.In(i);
        if strcmp(sg.kind, 'ref') && sg.dir > 3    % moment (world axes) on a reference point
            [~, Ra, ~, Jw] = s.anchorKin(s.Rf(sg.idx).anchor, q);
            e = zeros(3, 1);  e(sg.dir - 3) = u(i);
            F = F + Jw'*(Ra'*e);
        else                                       % joint force/torque, force on a point or node set
            [~, Jv] = s.coord(sg.kind, sg.idx, q);  F = F + Jv(sg.dir, :)'*u(i);
        end
    end
end
