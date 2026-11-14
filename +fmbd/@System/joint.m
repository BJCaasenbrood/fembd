function name = joint(s, A, B, type, opts)
%JOINT Kinematic joint between reference points A and B (rigid bodies,
%   couplings, or 'ground' for A). The released directions are joint
%   coordinates, zero in the configuration in which the bodies are
%   defined. Returns the joint name.
%
%   type      'fixed' (default), 'revolute', 'prismatic', 'cylindrical',
%             'spherical', 'planar', 'free' (bushing, use Stiffness)
%   Point     joint location (default: reference point of B)
%   Axis      joint axis = local z (default [0 0 1]); revolute and
%             prismatic act about/along it, planar moves normal to it
%   Frame     3x3 [ex ey ez] instead of Axis
%   Free      1x6 logical mask [ux uy uz rx ry rz] in the joint frame,
%             overrides type
%   Stiffness spring on the joint coordinates (scalar or per coordinate),
%             relaxed at zero
%   Damping   viscous damper on the joint coordinates
%   Position  initial joint coordinates (default 0)
%   Velocity  initial joint velocities (default 0)
%   Name      default 'joint<k>'
    arguments
        s
        A
        B
        type (1,:) char {mustBeMember(type, {'fixed', 'revolute', 'prismatic', ...
            'cylindrical', 'spherical', 'planar', 'free'})} = 'fixed'
        opts.Point double = []
        opts.Axis (1,3) double = [0 0 1]
        opts.Frame double = []
        opts.Free = []
        opts.Stiffness double = 0
        opts.Damping double = 0
        opts.Position double = []
        opts.Velocity double = []
        opts.Name (1,:) char = ''
    end
    iA = s.refIndex(A);  iB = s.refIndex(B);
    assert(iB > 0 && iA ~= iB, 'fmbd:joint', 'B must be a reference point other than A.');
    mask = struct('fixed', [0 0 0 0 0 0], 'revolute', [0 0 0 0 0 1], ...
        'prismatic', [0 0 1 0 0 0], 'cylindrical', [0 0 1 0 0 1], ...
        'spherical', [0 0 0 1 1 1], 'planar', [1 1 0 0 0 1], 'free', [1 1 1 1 1 1]);
    free = logical(mask.(type));
    if ~isempty(opts.Free), free = logical(opts.Free(:)'); end
    nf = nnz(free);
    R = opts.Frame;
    if isempty(R)
        ez = opts.Axis/norm(opts.Axis);
        h = [1 0 0];  if abs(ez(1)) > 0.9, h = [0 1 0]; end
        ex = h - dot(h, ez)*ez;  ex = ex/norm(ex);
        R = [ex; cross(ez, ex); ez]';
    end
    p = opts.Point;  if isempty(p), p = s.Refs(iB).x; end
    for f = {'Position', 'Velocity'}
        v = opts.(f{1});
        assert(isempty(v) || numel(v) == nf, 'fmbd:joint', '%s: %d value(s) expected.', f{1}, nf);
    end
    name = opts.Name;
    if isempty(name), name = sprintf('joint%d', numel(s.Joints) + 1); end
    assert(~ismember(name, s.names()), 'fmbd:joint', 'The name %s is already used.', name);
    s.Joints(end+1) = struct('Name', name, 'A', iA, 'B', iB, 'p', p(:)', 'R', R, 'free', free, ...
        'k', opts.Stiffness(:).*ones(nf, 1), 'damp', opts.Damping(:).*ones(nf, 1), ...
        'position', opts.Position(:), 'velocity', opts.Velocity(:));
    s.built = false;
end
