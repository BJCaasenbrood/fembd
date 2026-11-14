function setMass(b, m, opts)
%SETMASS Set the mass properties of a rigid body (e.g. from CAD).
%   b.setMass(m) sets the mass and scales the inertia (same distribution).
%   Options:
%     CenterOfMass  1x3 centre of mass (default: keep the current one)
%     Inertia       3x3 inertia tensor, or 1x3 principal moments
%                   (default: current inertia scaled by the mass ratio)
%     InertiaPoint  point about which Inertia is given (default: the
%                   centre of mass), e.g. the CAD part origin
%     Axes          3x3 axes in which Inertia is given (columns, in the
%                   assembly frame; default eye(3)), e.g. principal axes
%     LengthScale   factor for lengths in CenterOfMass, InertiaPoint and
%                   Inertia, e.g. 1e-3 for mm and kg mm^2 (default 1)
%   Units otherwise SI, in the assembly frame. A meshed rigid body uses
%   these values instead of its mesh mass.
    arguments
        b
        m (1,1) double {mustBePositive}
        opts.CenterOfMass double = []
        opts.Inertia double = []
        opts.InertiaPoint double = []
        opts.Axes (3,3) double = eye(3)
        opts.LengthScale (1,1) double {mustBePositive} = 1
    end
    assert(b.Rigid, 'fmbd:setMass', '%s: setMass is for rigid bodies only.', b.Name);
    p = struct();  [p.m, p.c, p.J] = b.massProperties();
    s = opts.LengthScale;  c = p.c;
    if ~isempty(opts.CenterOfMass), c = s*opts.CenterOfMass(:)'; end
    if isempty(opts.Inertia)
        J = p.J*(m/p.m);
    else
        J = opts.Inertia;
        if isvector(J), J = diag(J); end
        J = s^2*opts.Axes*J*opts.Axes';
        if ~isempty(opts.InertiaPoint) % parallel-axis theorem to the centre of mass
            r = c - s*opts.InertiaPoint(:)';
            J = J - m*(dot(r, r)*eye(3) - r'*r);
        end
    end
    assert(norm(J - J', 'fro') <= 1e-6*norm(J, 'fro'), 'fmbd:setMass', ...
        '%s: the inertia tensor must be symmetric.', b.Name);
    J = (J + J')/2;  Jp = sort(eig(J));
    assert(Jp(1) > 0 && Jp(3) <= (Jp(1) + Jp(2))*(1 + 1e-9), 'fmbd:setMass', ...
        ['%s: invalid inertia (principal moments %s): they must be positive and ' ...
         'satisfy J3 <= J1 + J2. Check that it is about the centre of mass.'], ...
        b.Name, mat2str(Jp', 4));
    b.MassProps = struct('m', m, 'c', c, 'J', J);
end
