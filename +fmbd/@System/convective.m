function gd = convective(s, q, nu)
% dG/dt * nu, by central differences along the motion.
    if norm(nu) == 0, [g, ~] = s.constraints(q);  gd = zeros(size(g)); return, end
    e = 1e-6/norm(nu);
    [~, Gp] = s.constraints(s.retract(q, e*nu));
    [~, Gm] = s.constraints(s.retract(q, -e*nu));
    gd = (Gp - Gm)*nu/(2*e);
end
