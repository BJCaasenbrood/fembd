function K = stiffness(s, q, u, lam, nu, cols)
% Tangent stiffness -d(F - G'*lam)/dq by central differences (columns cols).
    if nargin < 5 || isempty(nu), nu = zeros(s.nv, 1); end
    if nargin < 6, cols = 1:s.nv; end
    n = s.nv;
    K = zeros(n);
    e = 1e-5; % ~ eps^(1/3): truncation vs round-off
    for j = cols(:)'
        d = zeros(n, 1);  d(j) = e;
        qp = retract(s, q, d);  qm = retract(s, q, -d);
        [~, Gp] = constraints(s, qp);  [~, Gm] = constraints(s, qm);
        K(:, j) = -((forces(s, qp, nu, u, false) - Gp'*lam) - (forces(s, qm, nu, u, false) - Gm'*lam))/(2*e);
    end
end
