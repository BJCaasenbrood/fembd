function b = rigid(file, rho, opts)
%RIGID Rigid body from a closed STL surface.
%   b = fmbd.Body.rigid(file, rho) computes the mass, centre of mass and
%   inertia from the surface. Options: Name, Scale, RefPoint.
    arguments
        file (1,:) char
        rho (1,1) double
        opts.Name (1,:) char = ''
        opts.Scale (1,1) double = 1
        opts.RefPoint double = []
    end
    tr = stlread(file);
    X = opts.Scale*tr.Points;  F = tr.ConnectivityList;
    [m, c, J, flip] = surfaceMass(X, F, rho);
    if flip, F = F(:, [1 3 2]); end
    if isempty(opts.Name), [~, opts.Name] = fileparts(file); end
    if isempty(opts.RefPoint), opts.RefPoint = c; end
    b = fmbd.Body(struct('X', X, 'T', zeros(0, 4)), struct('rho', rho), ...
                 Name = opts.Name, Rigid = true, RefPoint = opts.RefPoint);
    b.F = F;
    b.MassProps = struct('m', m, 'c', c, 'J', J);
end
