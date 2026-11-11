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
%   Every body has a moving frame: a rigid body at its centre of mass; an
%   elastic body at its first coupling, with its deformation described by
%   Craig-Bampton modes in that frame (first coupling clamped, the other
%   couplings as interfaces): x = r + R*(X0 + Phi*eta). Joints and supports
%   are exact nonlinear constraints with Lagrange multipliers. Large
%   rotations are exact; the deformation is linear.
%
%   Model description (the system is built when it is first used, and
%   rebuilt after any change):
%     add, couple, joint, transmission, fix, input, output
%   Analysis:
%     simulate, energy, animate, plot
%
%   Conventions:
%     The configuration in which the bodies are defined is the relaxed
%     reference: joint coordinates, reference-point coordinates and
%     transmission stretches are zero there, and joint springs and elastic
%     transmissions are relaxed there. The initial configuration has the
%     joint coordinates set by Position (default 0); the other coordinates
%     follow from the joints, rigid transmissions and elastic transmissions
%     (unstretched, unless all their coordinates are set). Elastic bodies
%     start undeformed. Joint rotations are rotation vectors in the joint
%     frame (for a revolute joint: the angle). Node-set signals are the
%     average displacement of the nodes (world axes).

properties
    Gravity (1,3) double = [0 0 0]      % gravity vector, e.g. [0 0 -9.81]
end

properties (SetAccess = private)
    Bodies = {}                                                       % fmbd.Body handles
    Refs = struct('Name', {}, 'body', {}, 'nodes', {}, 'x', {})       % reference points
    Joints = struct('Name', {}, 'A', {}, 'B', {}, 'p', {}, 'R', {}, 'free', {}, 'k', {}, ...
                    'damp', {}, 'position', {}, 'velocity', {})
    Transmissions = struct('Name', {}, 'terms', {}, 'coef', {}, 'k', {}, 'damp', {})
    Fixes = struct('kind', {}, 'idx', {}, 'dirs', {})
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
    Bd                      % bodies: elastic, m, J, c0, Xl, ne, Phi, Q, s, CR, Kc, Dc, Mee, Psi, cols, ecols
    An                      % anchors (body frames, couplings): body, frame, x0, Psi
    Rf                      % reference points: anchor, s (offset in the anchor frame), x0
    Jn                      % joints and ground supports: A, B (anchors), a, b, QA, QB, free, k, c
    Lk                      % supported joint coordinates: joint, dir, value
    Tr                      % transmissions: terms (kind, idx, dir), coef, k, c
    Nd                      % node-set signals: body, u0, Phi, x0
    In                      % input signals: kind, idx, dir
    Out                     % output signals: kind, idx, dir, type
    N = 0                   % number of velocities
    q0                      % initial configuration: r (3 x nb), R (3 x 3 x nb), e
    nu0                     % initial velocities
    keep                    % independent constraint rows
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

    %% Model description
    function name = add(s, b)
    %ADD Add a body (fmbd.Body); returns its name. A rigid body is also a
    %   reference point with the body's name (for joints and signals).
        arguments
            s
            b (1,1) fmbd.Body
        end
        assert(~any(cellfun(@(c) c == b, s.Bodies)), 'fmbd:add', '%s is already added.', b.Name);
        assert(~ismember(b.Name, s.names()), 'fmbd:add', 'The name %s is already used.', b.Name);
        s.Bodies{end+1} = b;  name = b.Name;
        if b.Rigid
            s.Refs(end+1) = struct('Name', b.Name, 'body', numel(s.Bodies), 'nodes', [], 'x', b.Ref);
        end
        s.built = false;
    end

    function name = couple(s, body, sel, opts)
    %COUPLE Rigid coupling (RBE2): the nodes of an elastic body selected by
    %   sel (function of the node coordinates, or node indices) move with a
    %   new reference point, e.g. a bore for a joint or a bolt. Returns its
    %   name. The first coupling of a body carries the body's frame.
    %   Options: Name, Point (default: centroid of the selected nodes).
        arguments
            s
            body
            sel
            opts.Name (1,:) char = ''
            opts.Point double = []
        end
        [kind, i] = s.lookup(body, false);
        assert(strcmp(kind, 'body') && ~s.Bodies{i}.Rigid, 'fmbd:couple', ...
            'Couplings are for elastic bodies; a rigid body is already a reference point.');
        b = s.Bodies{i};  nodes = b.nodes(sel);
        assert(~isempty(nodes), 'fmbd:couple', 'No nodes selected.');
        x = opts.Point;  if isempty(x), x = mean(b.X(nodes, :), 1); end
        name = opts.Name;
        if isempty(name), name = sprintf('%s_ref%d', b.Name, numel(s.Refs) + 1); end
        assert(~ismember(name, s.names()), 'fmbd:couple', 'The name %s is already used.', name);
        s.Refs(end+1) = struct('Name', name, 'body', i, 'nodes', nodes, 'x', x(:)');
        s.built = false;
    end

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

    function name = transmission(s, ids, coeffs, opts)
    %TRANSMISSION Belt, gear or cable ratio between coordinates: the stretch
    %   delta = sum_i c_i s_i of joint coordinates or reference-point
    %   coordinates (Dirs 1:6 = [ux uy uz rx ry rz], absolute).
    %     Stiffness = k   elastic: a spring k on delta (with Damping)
    %     (none)          rigid: delta = 0
    %   Dirs gives the coordinate of each entry (default 1). Examples (ratio
    %   N = theta_motor/theta_driven):
    %     s.transmission({'motor', 'shoulder'}, [1/N -1], Stiffness = k)
    %     s.transmission({'motor', 'lower'}, [1/N -1], Dirs = {1, 6})   % absolute rz
        arguments
            s
            ids (1,:) cell
            coeffs (1,:) double
            opts.Dirs = []
            opts.Stiffness double = []
            opts.Damping (1,1) double = 0
            opts.Name (1,:) char = ''
        end
        nj = numel(ids);
        assert(numel(coeffs) == nj, 'fmbd:transmission', 'One coefficient per entry.');
        dirs = ones(1, nj);
        if iscell(opts.Dirs), dirs = [opts.Dirs{:}];
        elseif ~isempty(opts.Dirs), dirs = opts.Dirs(:)'; end
        assert(numel(dirs) == nj, 'fmbd:transmission', 'One entry of Dirs per coordinate.');
        terms = struct('kind', {}, 'idx', {}, 'dir', {});
        for i = 1:nj
            [kind, k] = s.lookup(ids{i}, true);
            assert(ismember(kind, {'joint', 'ref'}), 'fmbd:transmission', ...
                '%s is not a joint, rigid body or coupling.', string(ids{i}));
            nd = 6;  if strcmp(kind, 'joint'), nd = nnz(s.Joints(k).free); end
            assert(dirs(i) >= 1 && dirs(i) <= nd, 'fmbd:transmission', ...
                '%s has %d coordinate(s).', string(ids{i}), nd);
            terms(i) = struct('kind', kind, 'idx', k, 'dir', dirs(i));
        end
        assert(isempty(opts.Stiffness) || isscalar(opts.Stiffness), 'fmbd:transmission', ...
            'Stiffness must be a scalar (or empty for a rigid transmission).');
        assert(~isempty(opts.Stiffness) || opts.Damping == 0, 'fmbd:transmission', ...
            'Damping needs an elastic transmission (give Stiffness, possibly 0).');
        name = opts.Name;
        if isempty(name), name = sprintf('transmission%d', numel(s.Transmissions) + 1); end
        s.Transmissions(end+1) = struct('Name', name, 'terms', terms, 'coef', coeffs, ...
                                        'k', opts.Stiffness, 'damp', opts.Damping);
        s.built = false;
    end

    function fix(s, id, dirs)
    %FIX Support: s.fix(ref, dirs) holds a reference point (rigid body or
    %   coupling) in the directions dirs (default 1:6, [ux uy uz rx ry rz],
    %   world axes); s.fix(joint, dirs) holds joint coordinates at their
    %   Position (default all).
        if nargin < 3, dirs = []; end
        [kind, k] = s.lookup(id, true);
        assert(ismember(kind, {'joint', 'ref'}), 'fmbd:fix', ['%s is not a joint or reference ' ...
            'point; support an elastic body through a coupling.'], string(id));
        s.Fixes(end+1) = struct('kind', kind, 'idx', k, 'dirs', dirs);
        s.built = false;
    end

    function input(s, name, id, sel, dirs)
    %INPUT Named input, one signal per direction:
    %   s.input(name, joint)             force/torque on the joint coordinates
    %   s.input(name, ref, [], dirs)     force (1:3) or moment (4:6, world
    %                                    axes) on a rigid body or coupling
    %   s.input(name, body, sel, dirs)   total force on a node set, split
    %                                    equally over the nodes
        if nargin < 4, sel = []; end
        if nargin < 5, dirs = []; end
        s.Inputs = [s.Inputs, s.signalSpec(name, id, sel, dirs, 'force')];
        s.built = false;
    end

    function output(s, name, id, sel, dirs, opts)
    %OUTPUT Named output (sensor), one signal per direction: joint
    %   coordinates, reference-point coordinates or the average displacement
    %   of a node set, selected as in input. Type: 'displacement' (default)
    %   or 'velocity'.
        arguments
            s
            name (1,:) char
            id
            sel = []
            dirs = []
            opts.Type (1,:) char {mustBeMember(opts.Type, {'displacement', 'velocity'})} = 'displacement'
        end
        s.Outputs = [s.Outputs, s.signalSpec(name, id, sel, dirs, opts.Type)];
        s.built = false;
    end

    %% Analysis
    function [t, y, x] = simulate(s, tspan, opts)
    %SIMULATE Time integration with the generalised-alpha method.
    %   [t, y, x] = s.simulate([t0 t1]) starts in the initial configuration
    %   (joint Position and Velocity) and returns the times t, the outputs y
    %   (struct, one field per output signal) and the states x, a struct with
    %     r   3 x nb x nt     frame positions (centre of mass of rigid bodies)
    %     R   3 x 3 x nb x nt frame rotations
    %     e   ne x nt         modal coordinates of the elastic bodies
    %     nu  nDof x nt       velocities [v (world); omega (body); deta] per body
    %   Positions satisfy the constraints at every step (Newton tolerance),
    %   velocities are projected onto them (mass-orthogonal).
    %
    %   Options:
    %     Step   time step (default 1e-3)
    %     Input  function u = f(t, y) giving the input signals (column, in
    %            the order of s.InputName) from the outputs y (struct) of
    %            the previous step: a sampled controller (default: zero)
    %     Rho    spectral radius at infinity, 0..1 (default 0.8): 1 has no
    %            numerical damping, lower values damp high frequencies
    %     Tol    Newton tolerance on the accelerations, relative (default 1e-8)
        arguments
            s
            tspan (1,2) double
            opts.Step (1,1) double {mustBePositive} = 1e-3
            opts.Input = []
            opts.Rho (1,1) double {mustBeInRange(opts.Rho, 0, 1)} = 0.8
            opts.Tol (1,1) double {mustBePositive} = 1e-8
        end
        s.build();
        nt = max(1, round(diff(tspan)/opts.Step));  h = diff(tspan)/nt;
        t = tspan(1) + h*(0:nt)';
        r = opts.Rho;  am = (2*r - 1)/(r + 1);  af = r/(r + 1);
        gm = 0.5 + af - am;  be = 0.25*(gm + 0.5)^2;  ca = (1 - af)/(1 - am);
        nb = numel(s.Bd);  n = s.N;

        q = s.q0;  nu = s.nu0;  u = zeros(numel(s.In), 1);
        y0 = s.outputs(q, nu);
        if ~isempty(opts.Input), u = opts.Input(t(1), y0); end
        [F, ~, ~, Mq] = s.forces(q, nu, u);
        [~, G] = s.constraints(q);  nc = size(G, 1);
        sol = [Mq, G'; G, zeros(nc)]\[F; -s.convective(q, nu)];
        wd = sol(1:n);  lam = sol(n+1:end);  acc = wd;

        x = struct('r', zeros(3, nb, nt+1), 'R', zeros(3, 3, nb, nt+1), ...
                   'e', zeros(numel(q.e), nt+1), 'nu', zeros(n, nt+1));
        Y = zeros(numel(s.Out), nt+1);
        x.r(:,:,1) = q.r;  x.R(:,:,:,1) = q.R;  x.e(:,1) = q.e;  x.nu(:,1) = nu;
        Y(:,1) = cell2mat(struct2cell(y0));
        yk = y0;
        for k = 1:nt
            if ~isempty(opts.Input), u = opts.Input(t(k), yk); end
            w = wd;  prev = inf;
            for it = 1:30
                an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
                q1 = s.retract(q, h*(nu + h*(0.5 - be)*acc + h*be*an));
                nu1 = nu + h*(1 - gm)*acc + h*gm*an;
                [F, Kt, Ct, Mq] = s.forces(q1, nu1, u);
                [g, G] = s.constraints(q1);
                res = [Mq*w - F + G'*lam; g/(h^2*be*ca)];
                S = [Mq + h*gm*ca*Ct + h^2*be*ca*Kt, G'; G, zeros(nc)];
                dx = -S\res;
                w = w + dx(1:n);  lam = lam + dx(n+1:end);
                e = norm(dx(1:n))/(1 + norm(w));
                if e <= opts.Tol || (e >= prev && e < 1e3*opts.Tol), break, end   % converged, or at round-off
                prev = e;
            end
            if it == 30 && e > 1e3*opts.Tol
                warning('fmbd:simulate', 'Newton did not converge at t = %g s; try a smaller Step.', t(k+1));
            end
            an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
            q = s.retract(q, h*(nu + h*(0.5 - be)*acc + h*be*an));
            nu = nu + h*(1 - gm)*acc + h*gm*an;
            if nc > 0                            % project onto G nu = 0 (mass-orthogonal)
                [~, G] = s.constraints(q);  Mq = s.massMatrix(q);
                p = [Mq, G'; G, zeros(nc)]\[zeros(n, 1); -G*nu];  nu = nu + p(1:n);
            end
            acc = an;  wd = w;
            yk = s.outputs(q, nu);
            x.r(:,:,k+1) = q.r;  x.R(:,:,:,k+1) = q.R;  x.e(:,k+1) = q.e;  x.nu(:,k+1) = nu;
            Y(:,k+1) = cell2mat(struct2cell(yk));
        end
        y = struct();
        for i = 1:numel(s.Out)
            y.(matlab.lang.makeValidName(s.Outputs(i).Name)) = Y(i,:)';
        end
    end

    function E = energy(s, x)
    %ENERGY Energies along a simulation: E = s.energy(x) with x from
    %   simulate gives a struct of row vectors: kinetic, gravity (potential
    %   of gravity), strain (elastic bodies), springs (joint springs and
    %   elastic transmissions) and total. Dampers dissipate the total.
        s.build();
        nt = size(x.nu, 2);  z = zeros(1, nt);
        E = struct('kinetic', z, 'gravity', z, 'strain', z, 'springs', z, 'total', z);
        g = s.Gravity(:);
        for k = 1:nt
            q = struct('r', x.r(:,:,k), 'R', x.R(:,:,:,k), 'e', x.e(:,k));  nu = x.nu(:,k);
            E.kinetic(k) = 0.5*nu'*s.massMatrix(q)*nu;
            for i = 1:numel(s.Bd)
                B = s.Bd(i);  S = zeros(3, 1);
                if B.elastic
                    eta = q.e(B.ecols);  S = B.s*[1; eta];
                    E.strain(k) = E.strain(k) + 0.5*eta'*B.Kc*eta;
                end
                E.gravity(k) = E.gravity(k) - g'*(B.m*q.r(:,i) + q.R(:,:,i)*S);
            end
            for j = 1:numel(s.Jn)
                J = s.Jn(j);
                if any(J.k), c = s.jointKin(j, q);  E.springs(k) = E.springs(k) + 0.5*J.k'*c.^2; end
            end
            for T = s.Tr(arrayfun(@(T) ~isempty(T.k), s.Tr))
                E.springs(k) = E.springs(k) + 0.5*T.k*s.stretch(T, q)^2;
            end
        end
        E.total = E.kinetic + E.gravity + E.strain + E.springs;
    end

    function [Phi, X0] = modeShapes(s, body)
    %MODESHAPES Nodal mode shapes of an elastic body: [Phi, X0] =
    %   s.modeShapes(body) with Phi (3n x ne, interleaved [ux1 uy1 uz1 ...])
    %   in the body frame and X0 (n x 3) the nodes relative to the frame
    %   origin (the first coupling), so the nodes are at r + R*(X0 + Phi*eta)
    %   with r, R and eta from simulate (x.r, x.R, x.e).
        s.build();
        [kind, i] = s.lookup(body, false);
        assert(strcmp(kind, 'body') && s.Bd(i).elastic, 'fmbd:modeShapes', 'Not an elastic body.');
        Phi = s.Bd(i).Phi;  X0 = s.Bd(i).Xl;
    end

    function h = plot(s, x, k, opts)
    %PLOT Draw the system: s.plot() in its initial configuration, or
    %   s.plot(x, k) at step k of a simulation (default: the last). Option
    %   Scale magnifies the elastic deformation. Returns one graphics handle
    %   per body (patch, or marker for a point mass).
        arguments
            s
            x = []
            k = []
            opts.Scale (1,1) double = 1
        end
        s.build();
        if isempty(x), x = s.state(s.q0, s.nu0); end
        if isempty(k), k = size(x.nu, 2); end
        nb = numel(s.Bd);
        V = arrayfun(@(i) s.vertices(i, x, k, opts.Scale), 1:nb, 'UniformOutput', false);
        h = drawBodies(s.Bodies, V, num2cell(lines(nb), 2)');
        if nargout == 0, clear h, end
    end

    function animate(s, t, x, opts)
    %ANIMATE Animate a simulation: s.animate(t, x) with t, x from simulate.
    %   Options: Rate (default 1) = simulated seconds per wall-clock
    %   second, so 1 is real time and 0.1 slow motion; Scale (default 1)
    %   magnifies the elastic deformation; Title.
        arguments
            s
            t (:,1) double
            x struct
            opts.Rate (1,1) double {mustBePositive} = 1
            opts.Scale (1,1) double = 1
            opts.Title (1,:) char = ''
        end
        s.build();
        nb = numel(s.Bd);
        V = @(k) arrayfun(@(i) s.vertices(i, x, k, opts.Scale), 1:nb, 'UniformOutput', false);
        cla;  h = s.plot(x, 1, Scale = opts.Scale);
        X = reshape(permute(x.r, [1 3 2]), 3, [])';           % all frame positions
        span = max([arrayfun(@(B) max([vecnorm(B.Xl, 2, 2); 0]), s.Bd), 1e-3]);
        axis(reshape([min(X, [], 1) - span; max(X, [], 1) + span], 1, []));
        playback(h, (t(end) - t(1))/opts.Rate, @(c) frame(c));
        function frame(c)
            k = find(t <= t(1) + opts.Rate*c, 1, 'last');
            setVertices(h, V(k));
            title(sprintf('%s  t = %.3f s', opts.Title, t(k)));
        end
    end
end

methods (Access = private)
    %% Description helpers
    function n = names(s)
    % All names in use (bodies, reference points, joints).
        n = [cellfun(@(b) b.Name, s.Bodies, 'UniformOutput', false), {s.Refs.Name}, {s.Joints.Name}];
    end

    function [kind, k] = lookup(s, id, whole)
    % Resolve id to a 'joint', 'ref' or 'body' index. Joints and reference
    % points take precedence when no node selection is given (whole).
        if isa(id, 'fmbd.Body')
            k = find(cellfun(@(b) b == id, s.Bodies));
            assert(isscalar(k), 'fmbd:System', '%s is not part of the system.', id.Name);
            kind = 'body';
            if whole && id.Rigid, k = find([s.Refs.body] == k & cellfun(@isempty, {s.Refs.nodes}));  kind = 'ref'; end
            return
        end
        if whole
            k = find(strcmp({s.Joints.Name}, id));  kind = 'joint';
            if isscalar(k), return, end
            k = find(strcmp({s.Refs.Name}, id));    kind = 'ref';
            if isscalar(k), return, end
        end
        k = find(strcmp(cellfun(@(b) b.Name, s.Bodies, 'UniformOutput', false), id));  kind = 'body';
        assert(isscalar(k), 'fmbd:System', 'Unknown name %s.', string(id));
        if whole && s.Bodies{k}.Rigid, [kind, k] = s.lookup(s.Bodies{k}, true); end
    end

    function k = refIndex(s, id)
    % Reference point index of id; 0 for ground.
        k = 0;
        if (ischar(id) || isstring(id)) && strcmpi(id, 'ground'), return, end
        [kind, k] = s.lookup(id, true);
        assert(strcmp(kind, 'ref'), 'fmbd:joint', ...
            'Joint ends must be rigid bodies, couplings (s.couple) or ''ground''.');
    end

    function sig = signalSpec(s, name, id, sel, dirs, type)
    % One signal per direction: name, or name_dir if there are several.
        used = [s.InputName, s.OutputName];
        assert(~any(strcmp(used, name) | startsWith(used, [name '_'])), 'fmbd:io', ...
            'Signal name %s is already used.', name);
        [kind, k] = s.lookup(id, isempty(sel));
        nodes = [];
        if strcmp(kind, 'body')
            nodes = s.Bodies{k}.nodes(sel);  if isempty(sel), nodes = (1:s.Bodies{k}.nNode)'; end
            assert(~isempty(nodes), 'fmbd:io', '%s: no nodes selected.', name);
            if isempty(dirs), dirs = 1:3; end
        elseif isempty(dirs)
            dirs = 1:6;  if strcmp(kind, 'joint'), dirs = 1:nnz(s.Joints(k).free); end
        end
        sig = struct('Name', {}, 'kind', {}, 'idx', {}, 'dir', {}, 'nodes', {}, 'type', {});
        for d = dirs(:)'
            nm = name;  if numel(dirs) > 1, nm = sprintf('%s_%d', name, d); end
            sig(end+1) = struct('Name', nm, 'kind', kind, 'idx', k, 'dir', d, ...
                                'nodes', nodes, 'type', type);                         %#ok<AGROW>
        end
    end

    %% Build: frames, anchors, constraints, initial configuration
    function build(s)
        if s.built, return, end
        nb = numel(s.Bodies);  off = 0;  eoff = 0;
        s.An = struct('body', {}, 'frame', {}, 'x0', {}, 'Psi', {});
        s.Bd = [];
        for i = 1:nb
            b = s.Bodies{i};
            if b.Rigid
                [mb, c, J] = b.massProperties();
                B = struct('elastic', false, 'm', mb, 'J', J, 'c0', c(:), 'Xl', b.X - c(:)', ...
                           'ne', 0, 'Phi', zeros(3*b.nNode, 0), 'Q', [], 's', [], 'CR', [], ...
                           'Kc', [], 'Dc', [], 'Mee', [], 'Psi', zeros(6, 0, 0), 'cpl', []);
            else
                cpl = find([s.Refs.body] == i & ~cellfun(@isempty, {s.Refs.nodes}));
                assert(~isempty(cpl), 'fmbd:System', ['%s: an elastic body needs at least one ' ...
                       'coupling (s.couple); its first coupling carries the body frame.'], b.Name);
                B = elasticBody(b, s.Refs(cpl));  B.cpl = cpl;
            end
            B.cols = off + (1:6 + B.ne);  B.ecols = eoff + (1:B.ne);
            off = off + 6 + B.ne;  eoff = eoff + B.ne;
            s.Bd = [s.Bd, B];
            s.An(i) = struct('body', i, 'frame', true, 'x0', zeros(3, 1), 'Psi', zeros(6, B.ne));
        end
        s.N = off;
        q = struct('r', reshape([s.Bd.c0], 3, nb), 'R', repmat(eye(3), 1, 1, nb), 'e', zeros(eoff, 1));

        % reference points: rigid bodies -> body anchor; couplings -> own anchor
        s.Rf = struct('anchor', {}, 's', {}, 'x0', {});
        for k = 1:numel(s.Refs)
            i = s.Refs(k).body;  x = s.Refs(k).x(:);  B = s.Bd(i);
            if isempty(s.Refs(k).nodes)
                s.Rf(k) = struct('anchor', i, 's', x - B.c0, 'x0', x);
            elseif B.cpl(1) == k                 % the first coupling carries the body frame
                s.Rf(k) = struct('anchor', i, 's', zeros(3, 1), 'x0', x);
            else
                s.An(end+1) = struct('body', i, 'frame', false, 'x0', x - B.c0, ...
                                     'Psi', B.Psi(:,:,B.cpl == k));
                s.Rf(k) = struct('anchor', numel(s.An), 's', zeros(3, 1), 'x0', x);
            end
        end

        % joints (zero at the defining configuration), supports of reference
        % points as joints to ground, supports of joints as locks
        s.Jn = struct('A', {}, 'B', {}, 'a', {}, 'b', {}, 'QA', {}, 'QB', {}, 'free', {}, 'k', {}, 'c', {});
        for J = s.Joints
            s.addJoint(J.A, J.B, J.p(:), J.R, J.free, J.k, J.damp);
        end
        s.Lk = struct('joint', {}, 'dir', {}, 'value', {});
        for f = s.Fixes
            dirs = f.dirs;
            if strcmp(f.kind, 'ref')
                if isempty(dirs), dirs = 1:6; end
                s.addJoint(0, f.idx, s.Refs(f.idx).x(:), eye(3), ~ismember(1:6, dirs), ...
                           zeros(6 - numel(unique(dirs)), 1), 0);
            else
                J = s.Joints(f.idx);  if isempty(dirs), dirs = 1:nnz(J.free); end
                for d = dirs(:)'
                    v = 0;  if ~isempty(J.position), v = J.position(d); end
                    s.Lk(end+1) = struct('joint', f.idx, 'dir', d, 'value', v);
                end
            end
        end
        s.Tr = struct('terms', {}, 'coef', {}, 'k', {}, 'c', {});
        for T = s.Transmissions
            s.Tr(end+1) = struct('terms', T.terms, 'coef', T.coef, 'k', T.k, 'c', T.damp);
        end

        % signals
        s.Nd = struct('body', {}, 'u0', {}, 'Phi', {}, 'x0', {});
        s.In = s.compileSignals(s.Inputs);  s.Out = s.compileSignals(s.Outputs);

        % initial configuration: joint positions, transmissions unstretched,
        % elastic bodies undeformed; then the independent constraint rows
        s.keep = [];
        q = s.initialConfiguration(q);
        [~, G] = s.constraints(q, false);
        if ~isempty(G)
            [~, Rq, p] = qr(G', 'vector');
            r = sum(abs(diag(Rq)) > 1e-10*max(abs(diag(Rq))));
            s.keep = sort(p(1:r));
        end
        s.q0 = q;
        s.nu0 = s.initialVelocities(q);
        s.built = true;
    end

    function addJoint(s, iA, iB, p, R, free, k, c)
    % Joint between reference points iA (0 = ground) and iB at point p with
    % frame R (world, defining configuration; anchors are world-aligned).
        aA = 0;  a = p;
        if iA > 0, aA = s.Rf(iA).anchor;  a = p - s.Rf(iA).x0 + s.Rf(iA).s; end
        aB = s.Rf(iB).anchor;  b = p - s.Rf(iB).x0 + s.Rf(iB).s;
        assert(aA ~= aB, 'fmbd:joint', 'A joint connects a reference point to itself.');
        kk = zeros(6, 1);  cc = zeros(6, 1);
        kk(free) = k(:).*ones(nnz(free), 1);  cc(free) = c(:).*ones(nnz(free), 1);
        s.Jn(end+1) = struct('A', aA, 'B', aB, 'a', a, 'b', b, 'QA', R, 'QB', R, ...
                             'free', logical(free), 'k', kk, 'c', cc);
    end

    function sig = compileSignals(s, list)
    % Signals -> (kind, idx, dir, type); node sets become average-displacement signals.
        sig = struct('kind', {list.kind}, 'idx', {list.idx}, 'dir', {list.dir}, 'type', {list.type});
        for c = find(strcmp({list.kind}, 'body'))
            i = list(c).idx;  B = s.Bd(i);  nodes = list(c).nodes(:);
            d = 3*nodes' - [2; 1; 0];  P = zeros(3, B.ne);
            if B.ne > 0, P = reshape(mean(reshape(B.Phi(d(:), :), 3, numel(nodes), B.ne), 2), 3, B.ne); end
            s.Nd(end+1) = struct('body', i, 'u0', mean(B.Xl(nodes, :), 1)', 'Phi', P, ...
                                 'x0', mean(s.Bodies{i}.X(nodes, :), 1)');
            sig(c).kind = 'node';  sig(c).idx = numel(s.Nd);
        end
    end

    function q = initialConfiguration(s, q)
    % Satisfy the constraints with the joint positions, and elastic
    % transmissions unstretched unless all their coordinates are set; the
    % positions are ramped up for large angles. Elastic bodies stay undeformed.
        tgt = zeros(0, 3);                       % joint, coordinate, value
        for j = 1:numel(s.Joints)
            v = s.Joints(j).position;
            for d = 1:numel(v), tgt(end+1, :) = [j, d, v(d)]; end %#ok<AGROW>
        end
        given = @(T) all(arrayfun(@(t) strcmp(t.kind, 'joint') && ...
                     ~isempty(s.Joints(t.idx).position), T.terms));
        relax = s.Tr(arrayfun(@(T) ~isempty(T.k) && ~given(T), s.Tr));
        rigid = [s.Bd.cols];  ec = [s.Bd.ecols];
        if ~isempty(ec)                          % keep the modal coordinates at zero
            el = arrayfun(@(B) B.cols(7:end), s.Bd, 'UniformOutput', false);
            rigid = setdiff(rigid, [el{:}]);
        end
        ns = max(1, ceil(max([abs(tgt(:,3)); 0])/0.3));
        for step = 1:ns
            for it = 1:50
                [g, G] = s.constraints(q, false);
                for i = 1:size(tgt, 1)
                    [v, Jv] = s.coord('joint', tgt(i,1), q);
                    g = [g; v(tgt(i,2)) - step/ns*tgt(i,3)];  G = [G; Jv(tgt(i,2), :)];   %#ok<AGROW>
                end
                for T = relax
                    [d, Jd] = s.stretch(T, q);  g = [g; d];  G = [G; Jd];               %#ok<AGROW>
                end
                if isempty(g) || max(abs(g)) < 1e-12, break, end
                dq = zeros(s.N, 1);  dq(rigid) = -lsqminnorm(G(:, rigid), g);
                q = s.retract(q, dq);
            end
        end
        assert(isempty(g) || max(abs(g)) < 1e-9, 'fmbd:System', ['The initial configuration ' ...
            'cannot satisfy the joints, transmissions and joint positions (residual %.1e).'], max(abs(g)));
    end

    function nu = initialVelocities(s, q)
    % Velocities satisfying the constraints with the joint velocities;
    % elastic bodies start without deformation rates.
        [~, G] = s.constraints(q);  A = G;  rhs = zeros(size(G, 1), 1);
        for j = 1:numel(s.Joints)
            v = s.Joints(j).velocity;
            if isempty(v), continue, end
            [~, Jv] = s.coord('joint', j, q);
            A = [A; Jv];  rhs = [rhs; v(:)];                                       %#ok<AGROW>
        end
        if all(rhs == 0), nu = zeros(s.N, 1); return, end
        ec = [s.Bd.ecols];  Z = zeros(numel(ec), s.N);
        for i = 1:numel(s.Bd)
            c = s.Bd(i).cols;  Z(s.Bd(i).ecols, c(7:end)) = eye(s.Bd(i).ne);
        end
        A = [A; Z];  rhs = [rhs; zeros(numel(ec), 1)];
        nu = lsqminnorm(A, rhs);
        if norm(A*nu - rhs) > 1e-8*max(1, norm(rhs))
            warning('fmbd:System', 'Joint velocities are not consistent with the constraints; using the closest.');
        end
    end

    %% Kinematics, constraints and forces
    function x = state(~, q, nu)
    % A one-step state struct as returned by simulate.
        x = struct('r', q.r, 'R', q.R, 'e', q.e, 'nu', nu);
    end

    function V = vertices(s, i, x, k, scale)
    % Node positions of body i at step k (deformation magnified by scale);
    % the centre of mass for a point mass.
        B = s.Bd(i);  U = B.Xl;
        if isempty(U), V = x.r(:,i,k)';  return, end
        if B.elastic, U = U + scale*reshape(B.Phi*x.e(B.ecols, k), 3, [])'; end
        V = U*x.R(:,:,i,k)' + x.r(:,i,k)';
    end

    function [ra, Ra, Jv, Jw] = anchorKin(s, k, q)
    % Position ra and rotation Ra of anchor k (0 = ground) with the
    % Jacobians of its velocity (world) and angular velocity (anchor frame).
        Jv = zeros(3, s.N);  Jw = zeros(3, s.N);
        if k == 0, ra = zeros(3, 1);  Ra = eye(3);  return, end
        A = s.An(k);  B = s.Bd(A.body);  c = B.cols;
        r = q.r(:, A.body);  R = q.R(:,:,A.body);
        Jv(:, c(1:3)) = eye(3);
        if A.frame
            ra = r;  Ra = R;  Jw(:, c(4:6)) = eye(3);
            return
        end
        eta = q.e(B.ecols);  ut = A.x0 + A.Psi(1:3,:)*eta;  psi = A.Psi(4:6,:)*eta;  E = expso3(psi);
        ra = r + R*ut;  Ra = R*E;
        Jv(:, c(4:6)) = -R*hat(ut);  Jv(:, c(7:end)) = R*A.Psi(1:3,:);
        Jw(:, c(4:6)) = E';          Jw(:, c(7:end)) = dexp(psi)*A.Psi(4:6,:);
    end

    function [c, Jc, g, Gr] = jointKin(s, j, q)
    % Joint coordinates c = [d; phi] (translation in the A frame, rotation
    % vector of B relative to A) with Jacobian Jc, and the rotation
    % constraints g (dot products of frame axes) with Jacobian Gr.
        J = s.Jn(j);
        [rA, RA, JvA, JwA] = s.anchorKin(J.A, q);
        [rB, RB, JvB, JwB] = s.anchorKin(J.B, q);
        D = rB + RB*J.b - rA - RA*J.a;
        JD = (JvB - RB*hat(J.b)*JwB) - (JvA - RA*hat(J.a)*JwA);
        Rrel = J.QA'*RA'*RB*J.QB;
        phi = logso3(Rrel);
        c = [J.QA'*RA'*D; phi];
        Jc = [J.QA'*(hat(RA'*D)*JwA + RA'*JD); dexpInv(phi)*(J.QB'*JwB - Rrel'*J.QA'*JwA)];
        fr = J.free(4:6);
        switch nnz(fr)
            case 0, pairs = [2 3; 3 1; 1 2];
            case 1, k = find(fr);  i = setdiff(1:3, k);  pairs = [i(:), [k; k]];
            case 2, k = find(~fr);  pairs = circshift([1 2 3], -k);  pairs = pairs(1:2);
            otherwise, pairs = zeros(0, 2);
        end
        g = zeros(size(pairs, 1), 1);  Gr = zeros(size(pairs, 1), s.N);
        for p = 1:size(pairs, 1)
            ua = J.QA(:, pairs(p,1));  wb = J.QB(:, pairs(p,2));
            u = RA*ua;  w = RB*wb;
            g(p) = u'*w;
            Gr(p, :) = -cross(RB'*u, wb)'*JwB - cross(RA'*w, ua)'*JwA;
        end
    end

    function [v, Jv] = coord(s, kind, idx, q)
    % Coordinates and Jacobian of a joint (its free coordinates), a
    % reference point ([displacement; rotation vector]) or a node set
    % (average displacement), all zero in the defining configuration.
        switch kind
            case 'joint'
                [c, Jc] = s.jointKin(idx, q);  f = s.Jn(idx).free;
                v = c(f);  Jv = Jc(f, :);
            case 'ref'
                R = s.Rf(idx);
                [ra, Ra, Ja, Jw] = s.anchorKin(R.anchor, q);
                phi = logso3(Ra);
                v = [ra + Ra*R.s - R.x0; phi];
                Jv = [Ja - Ra*hat(R.s)*Jw; dexpInv(phi)*Jw];
            otherwise                            % node set
                Ns = s.Nd(idx);  B = s.Bd(Ns.body);  c = B.cols;
                Rb = q.R(:,:,Ns.body);  ub = Ns.u0 + Ns.Phi*q.e(B.ecols);
                v = q.r(:, Ns.body) + Rb*ub - Ns.x0;
                Jv = zeros(3, s.N);
                Jv(:, c(1:3)) = eye(3);  Jv(:, c(4:6)) = -Rb*hat(ub);  Jv(:, c(7:end)) = Rb*Ns.Phi;
        end
    end

    function [d, Jd] = stretch(s, T, q)
    % Transmission stretch delta = sum c_i s_i.
        d = 0;  Jd = zeros(1, s.N);
        for i = 1:numel(T.terms)
            [v, Jv] = s.coord(T.terms(i).kind, T.terms(i).idx, q);
            d = d + T.coef(i)*v(T.terms(i).dir);  Jd = Jd + T.coef(i)*Jv(T.terms(i).dir, :);
        end
    end

    function [g, G] = constraints(s, q, reduce)
    % Constraint values g(q) = 0 and their Jacobian G (rows s.keep).
        if nargin < 3, reduce = true; end
        g = zeros(0, 1);  G = zeros(0, s.N);
        for j = 1:numel(s.Jn)
            [c, Jc, gr, Gr] = s.jointKin(j, q);  lt = ~s.Jn(j).free(1:3);
            g = [g; c(lt); gr];  G = [G; Jc(lt, :); Gr];                      %#ok<AGROW>
        end
        for T = s.Tr(arrayfun(@(T) isempty(T.k), s.Tr))
            [d, Jd] = s.stretch(T, q);  g = [g; d];  G = [G; Jd];             %#ok<AGROW>
        end
        for L = s.Lk
            [v, Jv] = s.coord('joint', L.joint, q);
            g = [g; v(L.dir) - L.value];  G = [G; Jv(L.dir, :)];              %#ok<AGROW>
        end
        if reduce, g = g(s.keep);  G = G(s.keep, :); end
    end

    function Mq = massMatrix(s, q)
    % Mass matrix M(q) (constant for rigid bodies).
        Mq = zeros(s.N);
        for i = 1:numel(s.Bd)
            B = s.Bd(i);  c = B.cols;
            if B.elastic
                Mq(c, c) = elasticTerms(B, q.R(:,:,i), q.e(B.ecols));
            else
                Mq(c, c) = blkdiag(B.m*eye(3), B.J);
            end
        end
    end

    function [F, Kt, Ct, Mq] = forces(s, q, nu, u)
    % Generalised forces F (gravity, inertia, elastic, springs, dampers,
    % inputs), approximate tangents Kt = -dF/dq, Ct = -dF/dnu, and M(q).
        n = s.N;  F = zeros(n, 1);  Kt = zeros(n);  Ct = zeros(n);  Mq = zeros(n);
        g = s.Gravity(:);
        for i = 1:numel(s.Bd)
            B = s.Bd(i);  c = B.cols;  w = nu(c(4:6));  R = q.R(:,:,i);
            if B.elastic
                eta = q.e(B.ecols);
                [Mb, Fb] = elasticTerms(B, R, eta, nu(c), g);
                Mq(c, c) = Mb;  F(c) = Fb;
                Cb = zeros(numel(c), numel(c));  e = 1e-7*max(1, norm(nu(c)));   % velocity tangent
                for k = 1:numel(c)
                    dv = nu(c);  dv(k) = dv(k) + e;
                    [~, Fp] = elasticTerms(B, R, eta, dv, g);
                    Cb(:, k) = (Fb - Fp)/e;
                end
                Ct(c, c) = Cb;
                Kt(c(7:end), c(7:end)) = B.Kc;
            else
                Mq(c, c) = blkdiag(B.m*eye(3), B.J);
                F(c(1:3)) = B.m*g;
                F(c(4:6)) = -cross(w, B.J*w);
                Ct(c(4:6), c(4:6)) = hat(w)*B.J - hat(B.J*w);
            end
        end
        for j = 1:numel(s.Jn)
            J = s.Jn(j);  f = J.free & (J.k' ~= 0 | J.c' ~= 0);
            if ~any(f), continue, end
            [c, Jc] = s.jointKin(j, q);  c = c(f);  Jc = Jc(f, :);  k = J.k(f);  d = J.c(f);
            F = F - Jc'*(k.*c + d.*(Jc*nu));
            Kt = Kt + Jc'*(k.*Jc);  Ct = Ct + Jc'*(d.*Jc);
        end
        for T = s.Tr(arrayfun(@(T) ~isempty(T.k), s.Tr))
            [d, Jd] = s.stretch(T, q);
            F = F - Jd'*(T.k*d + T.c*(Jd*nu));
            Kt = Kt + T.k*(Jd'*Jd);  Ct = Ct + T.c*(Jd'*Jd);
        end
        for i = 1:numel(s.In)
            sg = s.In(i);
            if strcmp(sg.kind, 'ref') && sg.dir > 3    % moment (world axes) on a reference point
                [~, Ra, ~, Jw] = s.anchorKin(s.Rf(sg.idx).anchor, q);
                e = zeros(3, 1);  e(sg.dir - 3) = u(i);
                F = F + Jw'*(Ra'*e);
            else                                       % joint force/torque, force on a point or node set
                [~, Jv] = s.coord(sg.kind, sg.idx, q);  F = F + Jv(sg.dir, :)'*u(i);
            end
        end
    end

    function y = outputs(s, q, nu)
    % Output signals as a struct (fields in the order of OutputName).
        v = zeros(numel(s.Out), 1);
        for i = 1:numel(s.Out)
            sg = s.Out(i);  [c, Jc] = s.coord(sg.kind, sg.idx, q);
            if strcmp(sg.type, 'velocity'), v(i) = Jc(sg.dir, :)*nu; else, v(i) = c(sg.dir); end
        end
        y = cell2struct(num2cell(v), matlab.lang.makeValidName(s.OutputName(:)), 1);
        if isempty(s.Out), y = struct(); end
    end

    function gd = convective(s, q, nu)
    % dG/dt * nu, by central differences along the motion.
        if norm(nu) == 0, [g, ~] = s.constraints(q);  gd = zeros(size(g)); return, end
        e = 1e-6/norm(nu);
        [~, Gp] = s.constraints(s.retract(q, e*nu));
        [~, Gm] = s.constraints(s.retract(q, -e*nu));
        gd = (Gp - Gm)*nu/(2*e);
    end

    function q = retract(s, q, dq)
    % q o exp(dq): positions add (world), rotations update on the right
    % (body), modal coordinates add.
        for i = 1:numel(s.Bd)
            c = s.Bd(i).cols;
            q.r(:, i) = q.r(:, i) + dq(c(1:3));
            q.R(:,:,i) = q.R(:,:,i)*expso3(dq(c(4:6)));
            q.e(s.Bd(i).ecols) = q.e(s.Bd(i).ecols) + dq(c(7:end));
        end
    end
end
end

function B = elasticBody(b, cpl)
% Floating-frame data of an elastic body: frame at its first coupling,
% Craig-Bampton modes (first coupling clamped, the others as 6-dof
% interfaces; coupled nodes follow their coupling rigidly, RBE2), and the
% inertia invariants Q_ab = U_a' Ms U_b with U = [X0, Phi].
nn = b.nNode;  nc = numel(cpl);
owner = zeros(nn, 1);
for c = 1:nc
    if any(owner(cpl(c).nodes))
        error('fmbd:System', '%s: couplings %s and %s share nodes.', ...
              b.Name, cpl(max(owner(cpl(c).nodes))).Name, cpl(c).Name);
    end
    owner(cpl(c).nodes) = c;
end
free = find(~owner);  nf = numel(free);
rows = reshape(3*free' - [2; 1; 0], [], 1);  I = rows;  Jc = (1:3*nf)';  V = ones(3*nf, 1);
for c = 2:nc                                       % u = u_c + theta_c x (x - x_c)
    p = cpl(c).nodes(:);  d = b.X(p, :) - cpl(c).x;  base = 3*nf + 6*(c - 2);  np = numel(p);
    r = 3*p' - [2; 1; 0];
    I = [I; r(:); reshape(repmat(3*p' - [2; 1; 0], 3, 1), [], 1)];                   %#ok<AGROW>
    Jc = [Jc; repmat(base + (1:3)', np, 1); reshape(repmat(base + [4 5 6], 3, np), [], 1)]; %#ok<AGROW>
    H = zeros(3, 3, np);                            % -hat(d) per node
    H(1,2,:) = d(:,3);  H(1,3,:) = -d(:,2);  H(2,1,:) = -d(:,3);
    H(2,3,:) = d(:,1);  H(3,1,:) = d(:,2);   H(3,2,:) = -d(:,1);
    V = [V; ones(3*np, 1); H(:)];                                                      %#ok<AGROW>
end
T = sparse(I, Jc, V, 3*nn, 3*nf + 6*(nc - 1));
[K, M] = b.matrices();  D = b.Damping(1)*M + b.Damping(2)*K;
Kr = T'*K*T;  Mr = T'*M*T;  Dr = T'*D*T;
ib = 3*nf + (1:6*(nc - 1))';
[Kc, Mc, Vr] = fmbd.craigBampton((Kr + Kr')/2, (Mr + Mr')/2, ib, b.Modes);
ne = size(Vr, 2);
Phi = full(T*Vr);
Psi = zeros(6, ne, nc);
for c = 2:nc, Psi(:,:,c) = Vr(3*nf + 6*(c - 2) + (1:6), :); end
x0 = cpl(1).x(:);  Xl = b.X - x0';
Ms = M(1:3:end, 1:3:end);
Uf = reshape(cat(3, Xl, permute(reshape(Phi, 3, nn, ne), [2 1 3])), nn, []);   % [X0, Phi_1, ...] per node
MU = Ms*Uf;
Qf = full(Uf'*MU);                                     % 3(1+ne) x 3(1+ne)
Q = permute(reshape(Qf, 3, 1 + ne, 3, 1 + ne), [1 3 2 4]);   % Q(:,:,a,b) = U_a' Ms U_b
s = reshape(full(sum(MU, 1)), 3, 1 + ne);              % first moments of U_a
CR = zeros(3, 1 + ne, ne);                             % cross(U_a, Phi_j) moments
for aa = 1:1 + ne
    for j = 1:ne, CR(:, aa, j) = cr(Q(:,:,aa,j+1)); end
end
Dc = Vr'*Dr*Vr;
B = struct('elastic', true, 'm', full(sum(Ms(:))), 'J', [], 'c0', x0, 'Xl', Xl, 'ne', ne, ...
           'Phi', Phi, 'Q', Q, 's', s, 'CR', CR, 'Kc', Kc, 'Dc', (Dc + Dc')/2, 'Mee', Mc, ...
           'Psi', Psi, 'cpl', []);
end

function [Mb, Fb] = elasticTerms(B, R, eta, nu, g)
% Mass matrix and generalised forces (gravity, centrifugal and Coriolis,
% elastic and damping) of an elastic body in its floating frame.
% Velocities nu = [v (world); omega (body); deta].
ne = B.ne;  xi = [1; eta];  S = B.s*xi;  S1 = B.s(:, 2:end);
Quu = reshape(reshape(B.Q, 9, []) * kron(xi, xi), 3, 3);
J = trace(Quu)*eye(3) - Quu;  J = (J + J')/2;
Mwe = reshape(reshape(permute(B.CR, [1 3 2]), 3*ne, 1 + ne)*xi, 3, ne);
Mb = [B.m*eye(3), -R*hat(S), R*S1; (-R*hat(S))', J, Mwe; (R*S1)', Mwe', B.Mee];
Mb = (Mb + Mb')/2;
if nargout < 2, return, end
w = nu(4:6);  de = nu(7:end);  gb = R'*g;
Quz = reshape(reshape(B.Q(:,:,:,2:end), 9, []) * kron(de, xi), 3, 3);
W1 = -Quu*hat(w);  W2 = Quz;                     % u x (w x y) = w (u.y) - y (u.w)
hv = R*(cross(w, cross(w, S)) + 2*cross(w, S1*de));
hw = w*trace(W1) - W1'*w + 2*(w*trace(W2) - W2'*w);
Qju = reshape(reshape(B.Q(:,:,2:end,:), 9*ne, 1 + ne)*xi, 3, 3, ne);
Qjz = reshape(reshape(B.Q(:,:,2:end,2:end), 9*ne, ne)*de, 3, 3, ne);
he = zeros(ne, 1);
for j = 1:ne
    he(j) = w'*Qju(:,:,j)*w - (w'*w)*trace(Qju(:,:,j)) - 2*w'*cr(Qjz(:,:,j));
end
Fb = [B.m*g - hv; cross(S, gb) - hw; S1'*gb - he - B.Kc*eta - B.Dc*de];
end

function v = cr(Q)
% sum_pq Ms_pq y_p x z_q from Q = sum_pq Ms_pq y_p z_q'.
v = [Q(2,3) - Q(3,2); Q(3,1) - Q(1,3); Q(1,2) - Q(2,1)];
end
