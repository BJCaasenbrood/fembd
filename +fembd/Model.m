classdef Model < handle
%MODEL Nonlinear multibody model of an assembly (large rotations).
%   m = fembd.Model(a) or m = a.nonlinear() builds a nonlinear model from
%   the assembly a. Every body gets a moving frame (position and rotation
%   matrix); every joint and support becomes an exact nonlinear constraint,
%   enforced with Lagrange multipliers. Joint springs and dampers,
%   transmissions, inputs and outputs are taken over from the assembly.
%
%   Elastic bodies use a floating frame of reference: the frame follows
%   the body's first coupling, and the deformation in that frame is
%   described by Craig-Bampton modes (first coupling clamped, the other
%   couplings as interfaces): x = r + R*(X0 + Phi*eta). Large rotations
%   are exact; the deformation is linear (no stiffening by centrifugal or
%   gravity preload).
%
%   [t, y, x] = m.simulate(tspan, ...)   time integration, see simulate
%   m.animate(t, x)                      animation of a simulation
%   E = m.energy(x)                      kinetic and potential energies
%
%   Options:
%     Gravity  1x3 gravity vector (default [0 0 0]), e.g. [0 0 -9.81]
%     Modes    fixed-interface modes per elastic body (default 10); the
%              static modes of the interfaces are always included
%
%   Conventions:
%     Joint coordinates are measured from the configuration in which the
%     joint was created (as in Assembly.pose): joint springs are relaxed
%     there. Rotations of a joint are its rotation vector in the joint
%     frame (for a revolute joint: the angle). Transmissions and joint
%     supports are relaxed (elastic) or held (rigid) at the initial
%     configuration. Reference-point coordinates are the displacement of
%     the point and the rotation vector of its frame, absolute, from the
%     initial configuration. Node-set signals are the average
%     displacement of the nodes (world axes).
%
%   Limitations: elastic bodies need at least one coupling (their frame);
%   ties and supports on nodes are not supported; displacement and
%   velocity outputs.

properties (SetAccess = private)
    Assembly                % the source assembly
    Gravity (1,3) double
    Bodies                  % per body: name, elastic, m, J, c0, Xl, cols, ecols, and modal data
    Anchors                 % body frames and couplings: body, frame, x0, Psi (6 x ne)
    Joints                  % name, A, B (anchors, 0 = ground), a, b, QA, QB, free, k, c
    Trans                   % transmissions: terms (kind, idx, dir), coef, k, c, d0
    Locks                   % supported joint coordinates: kind, idx, dir, value
    Refs                    % reference points: anchor, s (offset in the anchor frame), x0
    Nodes                   % node-set signals: body, u0 (mean local position), Phi (3 x ne), x0
    In                      % input signals: kind, idx, dir
    Out                     % output signals: kind, idx, dir, type
    InputName  cell
    OutputName cell
    nDof                    % number of velocities (6 per body + modes)
    q0                      % initial configuration: r (3 x nb), R (3 x 3 x nb), e (modal coordinates)
    keep                    % independent constraint rows
end

methods
    function m = Model(a, opts)
        arguments
            a fembd.Assembly
            opts.Gravity (1,3) double = [0 0 0]
            opts.Modes (1,1) double {mustBeInteger, mustBeNonnegative} = 10
        end
        m.Assembly = a;  m.Gravity = opts.Gravity;
        assert(isempty(a.Links), 'fembd:Model', 'Ties are not supported by the nonlinear model yet.');

        % bodies: rigid (frame at the centre of mass) or elastic (frame at
        % the first coupling, Craig-Bampton modes); one anchor per body frame
        nb = numel(a.Bodies);  off = 0;  eoff = 0;
        m.Anchors = struct('body', {}, 'frame', {}, 'x0', {}, 'Psi', {});
        for i = 1:nb
            b = a.Bodies{i};
            if b.Rigid
                [mb, c, J] = b.massProperties();
                B = struct('name', b.Name, 'elastic', false, 'm', mb, 'J', J, 'c0', c(:), ...
                           'Xl', b.X - c(:)', 'ne', 0, ...
                           'Phi', zeros(3*b.nNode, 0), 'Q', [], 's', [], 'CR', [], 'Kc', [], ...
                           'Dc', [], 'Mee', [], 'couplings', {{}}, 'Psi', zeros(6, 0, 0));
            else
                B = elasticBody(a, i, opts.Modes);
            end
            B.cols = off + (1:6 + B.ne);  B.ecols = eoff + (1:B.ne);
            off = off + 6 + B.ne;  eoff = eoff + B.ne;
            m.Bodies = [m.Bodies, B];
            m.Anchors(i) = struct('body', i, 'frame', true, 'x0', zeros(3, 1), 'Psi', zeros(6, B.ne));
        end
        m.nDof = off;
        m.q0 = struct('r', reshape([m.Bodies.c0], 3, nb), 'R', repmat(eye(3), 1, 1, nb), ...
                      'e', zeros(eoff, 1));

        % reference points: rigid bodies -> body anchor; couplings -> own anchor
        m.Refs = struct('anchor', {}, 's', {}, 'x0', {});
        for k = 1:numel(a.Refs)
            i = a.refBody(k);  x = a.refPoint(k);  x = x(:);  B = m.Bodies(i);
            if isempty(a.Refs(k).nodes)
                m.Refs(k) = struct('anchor', i, 's', x - B.c0, 'x0', x);
            else
                c = find(strcmp(B.couplings, a.Refs(k).Name));
                if c == 1                        % the first coupling carries the body frame
                    m.Refs(k) = struct('anchor', i, 's', zeros(3, 1), 'x0', x);
                else
                    m.Anchors(end+1) = struct('body', i, 'frame', false, 'x0', x - B.c0, ...
                                              'Psi', B.Psi(:,:,c));
                    m.Refs(k) = struct('anchor', numel(m.Anchors), 's', zeros(3, 1), 'x0', x);
                end
            end
        end

        % joints, with frames on A and B from the joint value (pose)
        m.Joints = struct('name', {}, 'A', {}, 'B', {}, 'a', {}, 'b', {}, 'QA', {}, ...
                          'QB', {}, 'free', {}, 'k', {}, 'c', {});
        for J = a.Joints
            [Rj, tj] = jointMotion(J, J.value);
            m.addJoint(J.Name, J.A, J.B, J.p(:), J.R, J.p(:) + tj, Rj*J.R, J.free, J.k, J.damp);
        end

        % supports: of a reference point -> joint to ground; of a joint -> lock
        m.Locks = struct('kind', {}, 'idx', {}, 'dir', {}, 'value', {});
        for f = a.Fixes
            [kind, k] = a.lookup(f.id, isempty(f.nodes));
            dirs = f.dirs;
            switch kind
                case 'ref'
                    if isempty(dirs), dirs = 1:6; end
                    x = a.refPoint(k);
                    m.addJoint(sprintf('fix_%s', a.Refs(k).Name), 0, k, x(:), eye(3), x(:), eye(3), ...
                               ~ismember(1:6, dirs), zeros(6 - numel(unique(dirs)), 1), 0);
                case 'joint'
                    if isempty(dirs), dirs = 1:nnz(a.Joints(k).free); end
                    for d = dirs(:)'
                        m.Locks(end+1) = struct('kind', 'joint', 'idx', k, 'dir', d, 'value', 0);
                    end
                otherwise
                    error('fembd:Model', ['Supports on body nodes are not supported by the ' ...
                          'nonlinear model; support a coupling instead.']);
            end
        end
        for l = 1:numel(m.Locks)
            v = m.coord(m.Locks(l).kind, m.Locks(l).idx, m.q0);
            m.Locks(l).value = v(m.Locks(l).dir);
        end

        % transmissions: rest (elastic) or held (rigid) value at the initial configuration
        m.Trans = struct('terms', {}, 'coef', {}, 'k', {}, 'c', {}, 'd0', {});
        for t = a.Transmissions
            terms = struct('kind', t.kind, 'idx', num2cell(t.idx), 'dir', num2cell(t.dirs));
            m.Trans(end+1) = struct('terms', terms, 'coef', t.c, 'k', t.k, 'c', t.damp, 'd0', 0);
            m.Trans(end).d0 = m.stretch(m.Trans(end), m.q0);
        end

        % inputs and outputs
        m.Nodes = struct('body', {}, 'u0', {}, 'Phi', {}, 'x0', {});
        [m.In, m.InputName] = m.signals(a.Inputs);
        [m.Out, m.OutputName] = m.signals(a.Outputs);
        assert(all(ismember({m.Out.type}, {'displacement', 'velocity'})), 'fembd:Model', ...
            'Acceleration outputs are not supported by the nonlinear model yet.');

        % independent constraints (drop redundant rows), consistent start
        m.keep = [];
        [g, G] = m.constraints(m.q0, false);
        if ~isempty(G)
            [~, Rq, p] = qr(G', 'vector');
            r = sum(abs(diag(Rq)) > 1e-10*max(abs(diag(Rq))));
            m.keep = sort(p(1:r));
        end
        q = m.q0;
        for it = 1:50
            [g, G] = m.constraints(q);
            if isempty(g) || max(abs(g)) < 1e-12, break, end
            q = m.retract(q, -G'*((G*G')\g));
        end
        assert(isempty(g) || max(abs(g)) < 1e-9, 'fembd:Model', ...
            'The initial configuration does not satisfy the joint constraints.');
        m.q0 = q;
    end

    function [t, y, x] = simulate(m, tspan, opts)
    %SIMULATE Time integration with the generalised-alpha method.
    %   [t, y, x] = m.simulate([t0 t1]) starts at rest in the initial
    %   configuration and returns the times t, the outputs y (struct, one
    %   field per output signal) and the states x, a struct with
    %     r   3 x nb x nt   frame positions (centre of mass of rigid bodies)
    %     R   3 x 3 x nb x nt  frame rotations
    %     e   ne x nt       modal coordinates of the elastic bodies
    %     nu  nDof x nt     velocities [v (world); omega (body); deta] per body
    %   Positions satisfy the constraints at every step (Newton tolerance),
    %   velocities are projected onto them (mass-orthogonal).
    %
    %   Options:
    %     Step   time step (default 1e-3)
    %     Input  function u = f(t, y) giving the input signals (column, in
    %            the order of m.InputName) from the outputs y (struct) of
    %            the previous step: a sampled controller (default: zero)
    %     Rates  struct of initial joint velocities, e.g.
    %            struct('shoulder', [0 0 2]) (default: at rest)
    %     Rho    spectral radius at infinity, 0..1 (default 0.8): 1 has no
    %            numerical damping, lower values damp high frequencies
    %     Tol    Newton tolerance on the accelerations, relative (default 1e-8)
        arguments
            m
            tspan (1,2) double
            opts.Step (1,1) double {mustBePositive} = 1e-3
            opts.Input = []
            opts.Rates struct = struct()
            opts.Rho (1,1) double {mustBeInRange(opts.Rho, 0, 1)} = 0.8
            opts.Tol (1,1) double {mustBePositive} = 1e-8
        end
        nt = max(1, round(diff(tspan)/opts.Step));  h = diff(tspan)/nt;
        t = tspan(1) + h*(0:nt)';
        r = opts.Rho;  am = (2*r - 1)/(r + 1);  af = r/(r + 1);
        gm = 0.5 + af - am;  be = 0.25*(gm + 0.5)^2;  ca = (1 - af)/(1 - am);
        nb = numel(m.Bodies);  n = m.nDof;  nIn = numel(m.InputName);

        q = m.q0;  nu = zeros(n, 1);
        if ~isempty(fieldnames(opts.Rates)), nu = m.initialRates(q, opts.Rates); end
        u = zeros(nIn, 1);
        y0 = m.outputs(q, nu);
        if ~isempty(opts.Input), u = opts.Input(t(1), y0); end
        [F, ~, ~, Mq] = m.forces(q, nu, u);
        [~, G] = m.constraints(q);  nc = size(G, 1);
        sol = [Mq, G'; G, zeros(nc)]\[F; -m.convective(q, nu)];
        wd = sol(1:n);  lam = sol(n+1:end);  acc = wd;

        x = struct('r', zeros(3, nb, nt+1), 'R', zeros(3, 3, nb, nt+1), ...
                   'e', zeros(numel(q.e), nt+1), 'nu', zeros(n, nt+1));
        Y = zeros(numel(m.OutputName), nt+1);
        x.r(:,:,1) = q.r;  x.R(:,:,:,1) = q.R;  x.e(:,1) = q.e;  x.nu(:,1) = nu;
        Y(:,1) = cell2mat(struct2cell(y0));
        yk = y0;
        for k = 1:nt
            if ~isempty(opts.Input), u = opts.Input(t(k), yk); end
            w = wd;  prev = inf;
            for it = 1:30
                an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
                q1 = m.retract(q, h*(nu + h*(0.5 - be)*acc + h*be*an));
                nu1 = nu + h*(1 - gm)*acc + h*gm*an;
                [F, Kt, Ct, Mq] = m.forces(q1, nu1, u);
                [g, G] = m.constraints(q1);
                res = [Mq*w - F + G'*lam; g/(h^2*be*ca)];
                S = [Mq + h*gm*ca*Ct + h^2*be*ca*Kt, G'; G, zeros(nc)];
                dx = -S\res;
                w = w + dx(1:n);  lam = lam + dx(n+1:end);
                e = norm(dx(1:n))/(1 + norm(w));
                if e <= opts.Tol || (e >= prev && e < 1e3*opts.Tol), break, end   % converged, or at round-off
                prev = e;
            end
            if it == 30 && e > 1e3*opts.Tol
                warning('fembd:simulate', 'Newton did not converge at t = %g s; try a smaller Step.', t(k+1));
            end
            an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
            q = m.retract(q, h*(nu + h*(0.5 - be)*acc + h*be*an));
            nu = nu + h*(1 - gm)*acc + h*gm*an;
            if nc > 0                            % project onto G nu = 0 (mass-orthogonal)
                [~, G] = m.constraints(q);  Mq = m.massMatrix(q);
                p = [Mq, G'; G, zeros(nc)]\[zeros(n, 1); -G*nu];  nu = nu + p(1:n);
            end
            acc = an;  wd = w;
            yk = m.outputs(q, nu);
            x.r(:,:,k+1) = q.r;  x.R(:,:,:,k+1) = q.R;  x.e(:,k+1) = q.e;  x.nu(:,k+1) = nu;
            Y(:,k+1) = cell2mat(struct2cell(yk));
        end
        y = struct();
        for i = 1:numel(m.OutputName)
            y.(matlab.lang.makeValidName(m.OutputName{i})) = Y(i,:)';
        end
    end

    function E = energy(m, x)
    %ENERGY Energies along a simulation: E = m.energy(x) with x from
    %   simulate gives a struct of row vectors: kinetic, gravity (potential
    %   of gravity), strain (elastic bodies), springs (joint springs and
    %   elastic transmissions) and total. Dampers dissipate the total.
        nt = size(x.nu, 2);  z = zeros(1, nt);
        E = struct('kinetic', z, 'gravity', z, 'strain', z, 'springs', z, 'total', z);
        g = m.Gravity(:);
        for k = 1:nt
            q = struct('r', x.r(:,:,k), 'R', x.R(:,:,:,k), 'e', x.e(:,k));  nu = x.nu(:,k);
            E.kinetic(k) = 0.5*nu'*m.massMatrix(q)*nu;
            for i = 1:numel(m.Bodies)
                B = m.Bodies(i);  S = zeros(3, 1);
                if B.elastic
                    eta = q.e(B.ecols);  S = B.s*[1; eta];
                    E.strain(k) = E.strain(k) + 0.5*eta'*B.Kc*eta;
                end
                E.gravity(k) = E.gravity(k) - g'*(B.m*q.r(:,i) + q.R(:,:,i)*S);
            end
            for j = 1:numel(m.Joints)
                J = m.Joints(j);
                if any(J.k), s = m.jointKin(j, q);  E.springs(k) = E.springs(k) + 0.5*J.k'*s.^2; end
            end
            for T = m.Trans(arrayfun(@(T) ~isempty(T.k), m.Trans))
                E.springs(k) = E.springs(k) + 0.5*T.k*m.stretch(T, q)^2;
            end
        end
        E.total = E.kinetic + E.gravity + E.strain + E.springs;
    end

    function animate(m, t, x, opts)
    %ANIMATE Animate a simulation: m.animate(t, x) with t, x from simulate.
    %   Options: Rate (default 1) = simulated seconds per wall-clock
    %   second, so 1 is real time and 0.1 slow motion; Scale (default 1)
    %   magnifies the elastic deformation; Title.
        arguments
            m
            t (:,1) double
            x struct
            opts.Rate (1,1) double {mustBePositive} = 1
            opts.Scale (1,1) double = 1
            opts.Title (1,:) char = ''
        end
        nb = numel(m.Bodies);  V = @(k) arrayfun(@(i) m.vertices(i, x, k, opts.Scale), 1:nb, 'UniformOutput', false);
        cla;  h = drawBodies(m.Assembly.Bodies, V(1), num2cell(lines(nb), 2)');
        X = reshape(permute(x.r, [1 3 2]), 3, [])';           % all frame positions
        span = max([arrayfun(@(B) max([vecnorm(B.Xl, 2, 2); 0]), m.Bodies), 1e-3]);
        axis(reshape([min(X, [], 1) - span; max(X, [], 1) + span], 1, []));
        playback(h, (t(end) - t(1))/opts.Rate, @(s) frame(s));
        function frame(s)
            k = find(t <= t(1) + opts.Rate*s, 1, 'last');
            setVertices(h, V(k));
            title(sprintf('%s  t = %.3f s', opts.Title, t(k)));
        end
    end
end

methods (Access = private)
    function V = vertices(m, i, x, k, scale)
    % Node positions of body i at step k (deformation magnified by scale);
    % the centre of mass for a point mass.
        B = m.Bodies(i);  U = B.Xl;
        if isempty(U), V = x.r(:,i,k)';  return, end
        if B.elastic, U = U + scale*reshape(B.Phi*x.e(B.ecols, k), 3, [])'; end
        V = U*x.R(:,:,i,k)' + x.r(:,i,k)';
    end

    function addJoint(m, name, iA, iB, pA, RA, pB, RB, free, k, c)
    % Joint between reference points iA (0 = ground) and iB: frame (pA, RA)
    % on A and (pB, RB) on B, in world coordinates at the initial
    % configuration (anchor frames are world-aligned there).
        aA = 0;  a = pA;
        if iA > 0, aA = m.Refs(iA).anchor;  a = pA - m.Refs(iA).x0 + m.Refs(iA).s; end
        aB = m.Refs(iB).anchor;  b = pB - m.Refs(iB).x0 + m.Refs(iB).s;
        assert(aA ~= aB, 'fembd:Model', '%s connects a reference point to itself.', name);
        kk = zeros(6, 1);  cc = zeros(6, 1);
        kk(free) = k(:).*ones(nnz(free), 1);  cc(free) = c(:).*ones(nnz(free), 1);
        m.Joints(end+1) = struct('name', name, 'A', aA, 'B', aB, 'a', a, 'b', b, ...
            'QA', RA, 'QB', RB, 'free', logical(free), 'k', kk, 'c', cc);
    end

    function [ra, Ra, Jv, Jw] = anchorKin(m, k, q)
    % Position ra and rotation Ra of anchor k (0 = ground) with the
    % Jacobians of its velocity (world) and angular velocity (anchor frame).
        n = m.nDof;  Jv = zeros(3, n);  Jw = zeros(3, n);
        if k == 0, ra = zeros(3, 1);  Ra = eye(3);  return, end
        A = m.Anchors(k);  B = m.Bodies(A.body);  c = B.cols;
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

    function [s, Js, g, Gr] = jointKin(m, j, q)
    % Joint coordinates s = [d; phi] (translation in the A frame, rotation
    % vector of B relative to A) with Jacobian Js, and the rotation
    % constraints g (dot products of frame axes) with Jacobian Gr.
        J = m.Joints(j);
        [rA, RA, JvA, JwA] = m.anchorKin(J.A, q);
        [rB, RB, JvB, JwB] = m.anchorKin(J.B, q);
        D = rB + RB*J.b - rA - RA*J.a;
        JD = (JvB - RB*hat(J.b)*JwB) - (JvA - RA*hat(J.a)*JwA);
        Rrel = J.QA'*RA'*RB*J.QB;
        phi = logso3(Rrel);
        s = [J.QA'*RA'*D; phi];
        Js = [J.QA'*(hat(RA'*D)*JwA + RA'*JD); dexpInv(phi)*(J.QB'*JwB - Rrel'*J.QA'*JwA)];
        fr = J.free(4:6);
        switch nnz(fr)
            case 0, pairs = [2 3; 3 1; 1 2];
            case 1, k = find(fr);  i = setdiff(1:3, k);  pairs = [i(:), [k; k]];
            case 2, k = find(~fr);  pairs = circshift([1 2 3], -k);  pairs = pairs(1:2);
            otherwise, pairs = zeros(0, 2);
        end
        g = zeros(size(pairs, 1), 1);  Gr = zeros(size(pairs, 1), m.nDof);
        for p = 1:size(pairs, 1)
            ua = J.QA(:, pairs(p,1));  wb = J.QB(:, pairs(p,2));
            u = RA*ua;  w = RB*wb;
            g(p) = u'*w;
            Gr(p, :) = -cross(RB'*u, wb)'*JwB - cross(RA'*w, ua)'*JwA;
        end
    end

    function [v, Jv] = coord(m, kind, idx, q)
    % Coordinates and Jacobian of a joint (its free coordinates), a
    % reference point ([displacement; rotation vector], absolute) or a node
    % set (average displacement).
        switch kind
            case 'joint'
                [s, Js] = m.jointKin(idx, q);  f = m.Assembly.Joints(idx).free;
                v = s(f);  Jv = Js(f, :);
            case 'ref'
                R = m.Refs(idx);
                [ra, Ra, Ja, Jw] = m.anchorKin(R.anchor, q);
                phi = logso3(Ra);
                v = [ra + Ra*R.s - R.x0; phi];
                Jv = [Ja - Ra*hat(R.s)*Jw; dexpInv(phi)*Jw];
            otherwise                            % node set
                N = m.Nodes(idx);  B = m.Bodies(N.body);  c = B.cols;
                Rb = q.R(:,:,N.body);  ub = N.u0 + N.Phi*q.e(B.ecols);
                v = q.r(:, N.body) + Rb*ub - N.x0;
                Jv = zeros(3, m.nDof);
                Jv(:, c(1:3)) = eye(3);  Jv(:, c(4:6)) = -Rb*hat(ub);  Jv(:, c(7:end)) = Rb*N.Phi;
        end
    end

    function [d, Jd] = stretch(m, T, q)
    % Transmission stretch delta = sum c_i s_i (minus its initial value).
        d = -T.d0;  Jd = zeros(1, m.nDof);
        for i = 1:numel(T.terms)
            [v, Jv] = m.coord(T.terms(i).kind, T.terms(i).idx, q);
            d = d + T.coef(i)*v(T.terms(i).dir);  Jd = Jd + T.coef(i)*Jv(T.terms(i).dir, :);
        end
    end

    function [g, G] = constraints(m, q, reduce)
    % Constraint values g(q) = 0 and their Jacobian G (rows m.keep).
        if nargin < 3, reduce = true; end
        g = zeros(0, 1);  G = zeros(0, m.nDof);
        for j = 1:numel(m.Joints)
            [s, Js, gr, Gr] = m.jointKin(j, q);  lt = ~m.Joints(j).free(1:3);
            g = [g; s(lt); gr];  G = [G; Js(lt, :); Gr];                      %#ok<AGROW>
        end
        for T = m.Trans(arrayfun(@(T) isempty(T.k), m.Trans))
            [d, Jd] = m.stretch(T, q);  g = [g; d];  G = [G; Jd];             %#ok<AGROW>
        end
        for L = m.Locks
            [v, Jv] = m.coord(L.kind, L.idx, q);
            g = [g; v(L.dir) - L.value];  G = [G; Jv(L.dir, :)];              %#ok<AGROW>
        end
        if reduce, g = g(m.keep);  G = G(m.keep, :); end
    end

    function Mq = massMatrix(m, q)
    % Mass matrix M(q) (constant for rigid bodies).
        Mq = zeros(m.nDof);
        for i = 1:numel(m.Bodies)
            B = m.Bodies(i);  c = B.cols;
            if B.elastic
                Mq(c, c) = elasticTerms(B, q.R(:,:,i), q.e(B.ecols));
            else
                Mq(c, c) = blkdiag(B.m*eye(3), B.J);
            end
        end
    end

    function [F, Kt, Ct, Mq] = forces(m, q, nu, u)
    % Generalised forces F (gravity, inertia, elastic, springs, dampers,
    % inputs), approximate tangents Kt = -dF/dq, Ct = -dF/dnu, and M(q).
        n = m.nDof;  F = zeros(n, 1);  Kt = zeros(n);  Ct = zeros(n);  Mq = zeros(n);
        g = m.Gravity(:);
        for i = 1:numel(m.Bodies)
            B = m.Bodies(i);  c = B.cols;  w = nu(c(4:6));  R = q.R(:,:,i);
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
        for j = 1:numel(m.Joints)
            J = m.Joints(j);  f = J.free & (J.k' ~= 0 | J.c' ~= 0);
            if ~any(f), continue, end
            [s, Js] = m.jointKin(j, q);  s = s(f);  Js = Js(f, :);  k = J.k(f);  c = J.c(f);
            F = F - Js'*(k.*s + c.*(Js*nu));
            Kt = Kt + Js'*(k.*Js);  Ct = Ct + Js'*(c.*Js);
        end
        for T = m.Trans(arrayfun(@(T) ~isempty(T.k), m.Trans))
            [d, Jd] = m.stretch(T, q);
            F = F - Jd'*(T.k*d + T.c*(Jd*nu));
            Kt = Kt + T.k*(Jd'*Jd);  Ct = Ct + T.c*(Jd'*Jd);
        end
        for i = 1:numel(m.In)
            s = m.In(i);
            if strcmp(s.kind, 'ref') && s.dir > 3    % moment (world axes) on a reference point
                [~, Ra, ~, Jw] = m.anchorKin(m.Refs(s.idx).anchor, q);
                e = zeros(3, 1);  e(s.dir - 3) = u(i);
                F = F + Jw'*(Ra'*e);
            else                                     % joint force/torque, force on a point or node set
                [~, Jv] = m.coord(s.kind, s.idx, q);  F = F + Jv(s.dir, :)'*u(i);
            end
        end
    end

    function y = outputs(m, q, nu)
    % Output signals as a struct (fields in the order of m.OutputName).
        v = zeros(numel(m.Out), 1);
        for i = 1:numel(m.Out)
            s = m.Out(i);  [c, Jc] = m.coord(s.kind, s.idx, q);
            if strcmp(s.type, 'velocity'), v(i) = Jc(s.dir, :)*nu; else, v(i) = c(s.dir); end
        end
        y = cell2struct(num2cell(v), matlab.lang.makeValidName(m.OutputName(:)), 1);
        if isempty(m.Out), y = struct(); end
    end

    function [specs, names] = signals(m, list)
    % Input/output selections -> signals (kind, index, direction); node
    % sets become average-displacement signals.
        a = m.Assembly;  sig = a.signals(list);  names = {sig.Name};
        specs = struct('kind', {sig.kind}, 'idx', {sig.idx}, 'dir', {sig.dir}, 'type', {sig.type});
        for c = find(strcmp({sig.kind}, 'body'))
            s = sig(c);  i = s.idx;  B = m.Bodies(i);  nodes = s.nodes(:);
            d = 3*nodes' - [2; 1; 0];  P = zeros(3, B.ne);
            if B.ne > 0, P = reshape(mean(reshape(B.Phi(d(:), :), 3, numel(nodes), B.ne), 2), 3, B.ne); end
            m.Nodes(end+1) = struct('body', i, 'u0', mean(B.Xl(nodes, :), 1)', 'Phi', P, ...
                                    'x0', mean(a.Bodies{i}.X(nodes, :), 1)');
            specs(c).kind = 'node';  specs(c).idx = numel(m.Nodes);
        end
    end

    function nu = initialRates(m, q, rates)
    % Velocities satisfying the constraints with the given joint rates.
        [~, G] = m.constraints(q);  A = G;  bvec = zeros(size(G, 1), 1);
        for f = fieldnames(rates)'
            [kind, k] = m.Assembly.lookup(f{1}, true);
            assert(strcmp(kind, 'joint'), 'fembd:simulate', 'Rates: %s is not a joint.', f{1});
            [~, Jv] = m.coord('joint', k, q);  v = rates.(f{1});
            assert(numel(v) == size(Jv, 1), 'fembd:simulate', ...
                'Rates: %s has %d coordinate(s).', f{1}, size(Jv, 1));
            A = [A; Jv];  bvec = [bvec; v(:)];                               %#ok<AGROW>
        end
        e = [m.Bodies.ecols];
        Z = zeros(numel(e), m.nDof);                % start without elastic velocities
        for i = 1:numel(m.Bodies)
            c = m.Bodies(i).cols;  ec = m.Bodies(i).ecols;
            Z(ec, c(7:end)) = eye(numel(ec));
        end
        A = [A; Z];  bvec = [bvec; zeros(numel(e), 1)];
        nu = lsqminnorm(A, bvec);
        if norm(A*nu - bvec) > 1e-8*max(1, norm(bvec))
            warning('fembd:simulate', 'Rates are not consistent with the constraints; using the closest.');
        end
    end

    function gd = convective(m, q, nu)
    % dG/dt * nu, by central differences along the motion.
        if norm(nu) == 0, [g, ~] = m.constraints(q);  gd = zeros(size(g)); return, end
        e = 1e-6/norm(nu);
        [~, Gp] = m.constraints(m.retract(q, e*nu));
        [~, Gm] = m.constraints(m.retract(q, -e*nu));
        gd = (Gp - Gm)*nu/(2*e);
    end

    function q = retract(m, q, dq)
    % q o exp(dq): positions add (world), rotations update on the right
    % (body), modal coordinates add.
        for i = 1:numel(m.Bodies)
            c = m.Bodies(i).cols;
            q.r(:, i) = q.r(:, i) + dq(c(1:3));
            q.R(:,:,i) = q.R(:,:,i)*expso3(dq(c(4:6)));
            q.e(m.Bodies(i).ecols) = q.e(m.Bodies(i).ecols) + dq(c(7:end));
        end
    end
end
end

function B = elasticBody(a, i, nModes)
% Floating-frame data of elastic body i: frame at its first coupling,
% Craig-Bampton modes (first coupling clamped, the others as interfaces),
% and the inertia invariants Q_ab = U_a' Ms U_b with U = [X0, Phi].
b = a.Bodies{i};
cpl = find(arrayfun(@(r) r.body == b && ~isempty(r.nodes), a.Refs));
assert(~isempty(cpl), 'fembd:Model', ['%s: an elastic body needs at least one coupling ' ...
       '(a.couple) in the nonlinear model; its first coupling carries the body frame.'], b.Name);
ta = fembd.Assembly();  ta.add(b);  names = {a.Refs(cpl).Name};
for c = 1:numel(cpl)
    ta.couple(b, a.Refs(cpl(c)).nodes, Name = names{c}, Point = a.Refs(cpl(c)).x);
end
ta.fix(names{1}, [], 1:6);
[K, M, T, D] = ta.system();
ib = cell2mat(cellfun(@(nm) ta.qdof(ta.dof(nm)), names(2:end)', 'UniformOutput', false));
if isempty(ib), ib = zeros(0, 1); end
[Kc, Mc, V] = fembd.craigBampton(K, M, ib, nModes);
U = T*V;  ne = size(V, 2);  nn = b.nNode;
Phi = full(U(1:3*nn, :));
Psi = zeros(6, ne, numel(cpl));
for c = 1:numel(cpl), Psi(:,:,c) = full(U(ta.dof(names{c}), :)); end
x0 = a.refPoint(cpl(1));  Xl = b.X - x0(:)';
[~, Mb] = b.matrices(MassOnly = true);  Ms = Mb(1:3:end, 1:3:end);
Uf = reshape(cat(3, Xl, permute(reshape(Phi, 3, nn, ne), [2 1 3])), nn, []);   % [X0, Phi_1, ...] per node
MU = Ms*Uf;
Qf = full(Uf'*MU);                                     % 3(1+ne) x 3(1+ne)
Q = permute(reshape(Qf, 3, 1 + ne, 3, 1 + ne), [1 3 2 4]);   % Q(:,:,a,b) = U_a' Ms U_b
s = reshape(full(sum(MU, 1)), 3, 1 + ne);              % first moments of U_a
CR = zeros(3, 1 + ne, ne);                             % cross(U_a, Phi_j) moments
for aa = 1:1 + ne
    for j = 1:ne, CR(:, aa, j) = cr(Q(:,:,aa,j+1)); end
end
B = struct('name', b.Name, 'elastic', true, 'm', full(sum(Ms(:))), 'J', [], 'c0', x0(:), ...
           'Xl', Xl, 'ne', ne, 'Phi', Phi, 'Q', Q, 's', s, 'CR', CR, ...
           'Kc', Kc, 'Dc', (V'*D*V + (V'*D*V)')/2, 'Mee', Mc);
B.couplings = names;  B.Psi = Psi;
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
