classdef System < handle
%SYSTEM Flexible multibody system: rigid and elastic bodies, joints, inputs.
%   s = fmbd.System(Gravity = [0 0 -9.81]);
%   s.add(fmbd.Body.rigid('link.stl', 7850, Scale = 1e-3));
%   s.add(fmbd.Body.elastic('tool.stl', alu, Size = 3, Modes = 10));
%   s.joint('ground', 'link', 'revolute', Name = 'shoulder', Position = 0.5);
%   s.joint('link', s.couple('tool', @(X) X(:,1) < 0.201), 'fixed');
%   s.input('tau', 'shoulder');  s.output('th', 'shoulder');
%   [t, y, x] = s.simulate([0 1], Input = @(t, y) 20*(1 - y.th));
%   s.animate(t, x);
%
%   Bodies have a moving frame; rigid bodies use the COM frame, elastic
%   bodies use the first coupling frame. Joints and supports impose exact
%   nonlinear constraints with Lagrange multipliers. Large rotations are
%   exact; deformation is linear.
%
%   Elastic bodies tied together are reduced as one group; a fixed group is
%   grounded. Bodies are defined in a relaxed reference; elastic bodies start
%   undeformed.

properties
    Gravity (1,3) double = [0 0 0]      % gravity vector, e.g. [0 0 -9.81]
end

properties (SetAccess = private)
    Bodies = {}  % fmbd.Body handles
    Refs = struct('Name', {}, 'body', {}, 'nodes', {}, 'x', {})% reference points
    Joints = struct('Name', {}, 'A', {}, 'B', {}, 'p', {}, 'R', {}, 'free', {}, 'k', {}, ...
                    'damp', {}, 'position', {}, 'velocity', {}) 
    Transmissions = struct('Name', {}, 'terms', {}, 'coef', {}, 'k', {}, 'damp', {}) 
    Ties = struct('master', {}, 'slave', {}, 'snodes', {}, 'mnodes', {}, 'N', {}, 'n', {}, 'A', {}, 'k', {})
    Fixes = struct('kind', {}, 'idx', {}, 'dirs', {}, 'nodes', {})
    Inputs = struct('Name', {}, 'kind', {}, 'idx', {}, 'dir', {}, 'nodes', {}, 'type', {})
    Outputs = struct('Name', {}, 'kind', {}, 'idx', {}, 'dir', {}, 'nodes', {}, 'type', {})
end

properties (Dependent)
    InputName               % input signal names, in the order of simulate's Input
    OutputName              % output signal names
    nDof                    % number of velocities (6 per body + elastic modes)
end

properties (Access = private)
    built = false
    Bd      % units (rigid bodies, elastic groups): elastic, m, J, c0, Xl, ne, Phi,
            % Q, s, CR, Kc, Dc, Mee, Psi, cpl, cols, ecols
    Ub      % unit of each body
    Uo      % node offset of each body in its unit
    An      % anchors (body frames, couplings): body, frame, x0, Psi
    Rf      % reference points: anchor, s (offset in the anchor frame), x0
    Jn      % joints and ground supports: A, B (anchors), a, b, QA, QB, free, k, c
    Lk      % supported joint coordinates: joint, dir, value
    Tr      % transmissions: terms (kind, idx, dir), coef, k, c
    Nd      % node-set signals: body, u0, Phi, x0
    In      % input signals: kind, idx, dir
    Out     % output signals: kind, idx, dir, type
    N = 0   % number of velocities
    q0      % initial configuration: r (3 x nb), R (3 x 3 x nb), e
    nu0     % initial velocities
    keep    % independent constraint rows
end

methods
    function s = System(opts)
        arguments
            opts.Gravity (1,3) double = [0 0 0]
        end
        s.Gravity = opts.Gravity;
    end

    function names = get.InputName(s),  names = {s.Inputs.Name}; end
    function names = get.OutputName(s), names = {s.Outputs.Name}; end
    function n = get.nDof(s),  s.build();  n = s.N; end
end

methods
    name = add(s, b)
    name = couple(s, body, sel, opts)
    ns = tie(s, master, slave, opts)
    name = joint(s, A, B, type, opts)
    name = transmission(s, ids, coeffs, opts)
    fix(s, id, varargin)
    input(s, name, id, sel, dirs)
    output(s, name, id, sel, dirs, opts)
    update(s)
    [t, y, x] = simulate(s, tspan, opts)
    E = energy(s, x)
    [x, info] = equilibrium(s, opts)
    lin = linearize(s, x, opts)
    sys = statespace(s, x, opts)
    [f, V] = modes(s, x, n, opts)
    animateMode(s, x, v, opts)
    h = plotMode(s, x, v, opts)
    x = perturb(s, x, dq)
    [Phi, X0] = modeShapes(s, body)
    h = plot(s, x, k, opts)
    animate(s, t, x, opts)
end

methods (Access = private)
    n = names(s)
    [kind, k] = lookup(s, id, whole)
    k = refIndex(s, id)
    sig = signalSpec(s, name, id, sel, dirs, type)
    build(s)
    addJoint(s, iA, iB, p, R, free, k, c)
    sig = compileSignals(s, list)
    q = initialConfiguration(s, q)
    nu = initialVelocities(s, q)
    [q, nu] = startState(s, x)
    u = inputVector(s, u)
    K = stiffness(s, q, u, lam, nu, cols)
    K = frameStiffness(s, q, u, lam)
    [x, a, V0] = modeScale(s, x, v, amplitude)
    x = state(s, q, nu, u)
    V = vertices(s, i, x, k, scale)
    B = elasticGroup(s, members)
    [ra, Ra, Jv, Jw] = anchorKin(s, k, q)
    [c, Jc, g, Gr] = jointKin(s, j, q)
    [v, Jv] = coord(s, kind, idx, q)
    [d, Jd] = stretch(s, T, q)
    [g, G] = constraints(s, q, reduce)
    Mq = massMatrix(s, q)
    [F, Kt, Ct, Mq] = forces(s, q, nu, u, tangents)
    y = outputs(s, q, nu, a)
    gd = convective(s, q, nu)
    q = retract(s, q, dq)
end
end
