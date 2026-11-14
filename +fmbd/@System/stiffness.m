function K = stiffness(s, q, u, lam, nu, cols)
% Tangent stiffness K = -d(F - G'*lam)/dq by central differences on the
% tangent space (gravity, springs, elastic and constraint-force terms),
% for the columns cols (default all).
    if nargin < 5 || isempty(nu), nu = zeros(s.nv, 1); end
    if nargin < 6, cols = 1:s.nv; end
    n = s.nv;
    K = zeros(n);
    e = 1e-5; % ~ eps^(1/3): truncation vs round-off
    for j = cols(:)'
        d = zeros(n, 1);  d(j) = e;
        qp = s.retract(q, d);  qm = s.retract(q, -d);
        [~, Gp] = s.constraints(qp);  [~, Gm] = s.constraints(qm);
        K(:, j) = -((s.forces(qp, nu, u, false) - Gp'*lam) - (s.forces(qm, nu, u, false) - Gm'*lam))/(2*e);
    end
end
