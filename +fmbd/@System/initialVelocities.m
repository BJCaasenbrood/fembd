function nu = initialVelocities(s, q)
% Velocities satisfying the constraints with the joint velocities;
% elastic bodies start without deformation rates.
    [~, G] = s.constraints(q);
    A = G;
    rhs = zeros(size(G, 1), 1);
    for j = 1:numel(s.Joints)
        v = s.Joints(j).velocity;
        if isempty(v), continue, end
        [~, Jv] = s.coord('joint', j, q);
        A = [A; Jv];  rhs = [rhs; v(:)]; %#ok<AGROW>
    end
    if all(rhs == 0), nu = zeros(s.nv, 1); return, end
    ec = [s.units.ecols];  Z = zeros(numel(ec), s.nv);
    for i = 1:numel(s.units)
        c = s.units(i).cols;  Z(s.units(i).ecols, c(7:end)) = eye(s.units(i).ne);
    end
    A = [A; Z];  rhs = [rhs; zeros(numel(ec), 1)];
    nu = lsqminnorm(A, rhs);
    if norm(A*nu - rhs) > 1e-8*max(1, norm(rhs))
        warning('fmbd:System', 'Joint velocities are not consistent with the constraints; using the closest.');
    end
end
