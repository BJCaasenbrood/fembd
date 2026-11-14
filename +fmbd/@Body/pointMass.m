function b = pointMass(m, c, J, opts)
%POINTMASS Rigid body given by its mass properties only: mass m,
%   centre of mass c (1x3) and inertia tensor J (3x3, about c; default
%   0). It has no surface; connect it with a joint (System.joint).
%   Options: Name (default 'mass'), RefPoint (default c).
    arguments
        m (1,1) double
        c (1,3) double
        J (3,3) double = zeros(3)
        opts.Name (1,:) char = 'mass'
        opts.RefPoint double = []
    end
    if isempty(opts.RefPoint), opts.RefPoint = c; end
    b = fmbd.Body(struct('X', zeros(0, 3), 'T', zeros(0, 4)), struct('rho', 0), ...
                 Name = opts.Name, Rigid = true, RefPoint = opts.RefPoint);
    b.MassProps = struct('m', m, 'c', c, 'J', (J + J')/2);
end
