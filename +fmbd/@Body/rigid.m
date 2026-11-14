function b = rigid(file, rho, opts)
%RIGID Rigid body from a closed STL surface, without volume mesh.
%   b = fmbd.Body.rigid(file, rho) reads the STL triangles, which are
%   used for plotting; mass, centre of mass and inertia are computed
%   exactly from the surface (divergence theorem). Use setMass to
%   override them (e.g. with CAD values).
%   Options: Name, Scale, RefPoint (default: centre of mass).
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
