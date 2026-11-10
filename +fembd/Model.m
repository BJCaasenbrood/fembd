classdef Model < handle
%MODEL Nonlinear multibody model of an assembly (large rotations).
%   m = fembd.Model(a) or m = a.nonlinear() builds a nonlinear model from
%   the assembly a. Every rigid body gets a moving frame (position of its
%   centre of mass and a rotation matrix); every joint and support becomes
%   an exact nonlinear constraint, enforced with Lagrange multipliers.
%   Joint springs and dampers, transmissions, inputs and outputs are taken
%   over from the assembly.
%
%   [t, y, x] = m.simulate(tspan, ...)   time integration, see simulate
%   m.animate(t, x)                      animation of a simulation
%
%   Options:
%     Gravity  1x3 gravity vector (default [0 0 0]), e.g. [0 0 -9.81]
%
%   Conventions:
%     Joint coordinates are measured from the configuration in which the
%     joint was created (as in Assembly.pose): joint springs are relaxed
%     there. Rotations of a joint are its rotation vector in the joint
%     frame (for a revolute joint: the angle). Transmissions and joint
%     supports are relaxed (elastic) or held (rigid) at the initial
%     configuration. Reference-point coordinates are the displacement of
%     the point and the rotation vector of the body, absolute, from the
%     initial configuration.
%
%   Limitations: rigid bodies only (elastic bodies, couplings and ties are
%   not supported yet); displacement and velocity outputs.

properties (SetAccess = private)
    Assembly                % the source assembly
    Gravity (1,3) double
    Bodies                  % name, m, J (about c, body frame), c0, Xl (local nodes), F
    Joints                  % name, A, B (body indices, 0 = ground), a, b, QA, QB, free, k, c
    Trans                   % transmissions: terms (kind, idx, dir), coef, k, c, d0
    Locks                   % supported joint coordinates: kind, idx, dir, value
    Refs                    % reference points: body, s (local offset), x0
    In                      % input signals: kind, idx, dir
    Out                     % output signals: kind, idx, dir, type
    InputName  cell
    OutputName cell
    M                       % constant mass matrix (6 per body: [m I, J])
    q0                      % initial configuration: r (3 x nb), R (3 x 3 x nb)
    keep                    % independent constraint rows
end

methods
    function m = Model(a, opts)
        arguments
            a fembd.Assembly
            opts.Gravity (1,3) double = [0 0 0]
        end
        m.Assembly = a;  m.Gravity = opts.Gravity;
        assert(all(cellfun(@(b) b.Rigid, a.Bodies)), 'fembd:Model', ...
            'Elastic bodies are not supported by the nonlinear model yet.');
        assert(isempty(a.Links), 'fembd:Model', 'Ties are not supported by the nonlinear model yet.');
        assert(all(arrayfun(@(r) isempty(r.nodes), a.Refs)), 'fembd:Model', ...
            'Couplings are not supported by the nonlinear model yet.');

        % bodies: mass properties, frame at the centre of mass (world-aligned)
        nb = numel(a.Bodies);
        m.Bodies = struct('name', {}, 'm', {}, 'J', {}, 'c0', {}, 'Xl', {}, 'F', {});
        for i = 1:nb
            b = a.Bodies{i};
            [mb, c, J] = massProps(b);
            m.Bodies(i) = struct('name', b.Name, 'm', mb, 'J', J, 'c0', c(:), ...
                                 'Xl', b.X - c(:)', 'F', b.F(:, 1:min(3, size(b.F, 2))));
        end
        m.M = zeros(6*nb);
        for i = 1:nb
            m.M(6*i-5:6*i, 6*i-5:6*i) = blkdiag(m.Bodies(i).m*eye(3), m.Bodies(i).J);
        end
        m.q0 = struct('r', [m.Bodies.c0], 'R', repmat(eye(3), 1, 1, nb));
        if nb == 0, m.q0.r = zeros(3, 0); end

        % reference points
        m.Refs = struct('body', {}, 's', {}, 'x0', {});
        for k = 1:numel(a.Refs)
            i = a.refBody(k);  x = a.refPoint(k);
            m.Refs(k) = struct('body', i, 's', x(:) - m.Bodies(i).c0, 'x0', x(:));
        end

        % joints, with frames on A and B from the joint value (pose)
        m.Joints = struct('name', {}, 'A', {}, 'B', {}, 'a', {}, 'b', {}, 'QA', {}, ...
                          'QB', {}, 'free', {}, 'k', {}, 'c', {});
        for J = a.Joints
            s = zeros(6, 1);  s(J.free) = J.value;
            m.addJoint(J.Name, J.A, J.B, J.p(:), J.R, J.p(:) + J.R*s(1:3), J.R*expso3(s(4:6)), ...
                       J.free, J.k, J.damp);
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
                    error('fembd:Model', 'Supports on body nodes are not supported by the nonlinear model.');
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
    %   field per output signal) and the states x (struct with r: 3 x nb x
    %   nt centre-of-mass positions, R: 3 x 3 x nb x nt rotations and nu:
    %   6nb x nt velocities [v (world); omega (body)] per body).
    %
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
        nb = numel(m.Bodies);  n = 6*nb;  nu0 = zeros(n, 1);
        Mm = m.M;  nIn = numel(m.InputName);

        q = m.q0;
        if ~isempty(fieldnames(opts.Rates)), nu0 = m.initialRates(q, opts.Rates); end
        nu = nu0;
        u = zeros(nIn, 1);
        y0 = m.outputs(q, nu);
        if ~isempty(opts.Input), u = opts.Input(t(1), y0); end
        [F, ~, ~] = m.forces(q, nu, u);
        [~, G] = m.constraints(q);  nc = size(G, 1);
        S = [Mm, G'; G, zeros(nc)];
        sol = S\[F; -m.convective(q, nu)];
        wd = sol(1:n);  lam = sol(n+1:end);  acc = wd;

        x = struct('r', zeros(3, nb, nt+1), 'R', zeros(3, 3, nb, nt+1), 'nu', zeros(n, nt+1));
        Y = zeros(numel(m.OutputName), nt+1);
        x.r(:,:,1) = q.r;  x.R(:,:,:,1) = q.R;  x.nu(:,1) = nu;  Y(:,1) = cell2mat(struct2cell(y0));
        yk = y0;
        for k = 1:nt
            if ~isempty(opts.Input), u = opts.Input(t(k), yk); end
            w = wd;  prev = inf;
            for it = 1:30
                an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
                q1 = m.retract(q, h*(nu + h*(0.5 - be)*acc + h*be*an));
                nu1 = nu + h*(1 - gm)*acc + h*gm*an;
                [F, Kt, Ct] = m.forces(q1, nu1, u);
                [g, G] = m.constraints(q1);
                res = [Mm*w - F + G'*lam; g/(h^2*be*ca)];
                S = [Mm + h*gm*ca*Ct + h^2*be*ca*Kt, G'; G, zeros(nc)];
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
                [~, G] = m.constraints(q);
                p = [Mm, G'; G, zeros(nc)]\[zeros(n, 1); -G*nu];  nu = nu + p(1:n);
            end
            acc = an;  wd = w;
            yk = m.outputs(q, nu);
            x.r(:,:,k+1) = q.r;  x.R(:,:,:,k+1) = q.R;  x.nu(:,k+1) = nu;
            Y(:,k+1) = cell2mat(struct2cell(yk));
        end
        y = struct();
        for i = 1:numel(m.OutputName)
            y.(matlab.lang.makeValidName(m.OutputName{i})) = Y(i,:)';
        end
    end

    function animate(m, t, x, opts)
    %ANIMATE Animate a simulation: m.animate(t, x) with t, x from simulate.
    %   Options: Rate (default 1) = simulated seconds per wall-clock
    %   second, so 1 is real time and 0.1 slow motion; Title.
        arguments
            m
            t (:,1) double
            x struct
            opts.Rate (1,1) double {mustBePositive} = 1
            opts.Title (1,:) char = ''
        end
        nb = numel(m.Bodies);  h = gobjects(nb, 1);
        hold on
        for i = 1:nb
            B = m.Bodies(i);
            if isempty(B.F)                      % point mass: marker at its centre of mass
                h(i) = plot3(x.r(1,i,1), x.r(2,i,1), x.r(3,i,1), 'ko', 'MarkerFaceColor', 'k');
            else
                h(i) = patch('Faces', B.F, 'Vertices', B.Xl*x.R(:,:,i,1)' + x.r(:,i,1)', ...
                             'FaceColor', [0.75 0.78 0.85], 'LineStyle', 'none');
            end
        end
        hold off;  axis equal;  view(3);  camlight;  lighting gouraud
        X = reshape(permute(x.r, [1 3 2]), 3, [])';           % all centre-of-mass positions
        span = max([arrayfun(@(B) max([vecnorm(B.Xl, 2, 2); 0]), m.Bodies), 1e-3]);
        axis(reshape([min(X, [], 1) - span; max(X, [], 1) + span], 1, []));
        clock = tic;  k = 1;
        while k < numel(t) && all(isgraphics(h))
            k = find(t <= t(1) + opts.Rate*toc(clock), 1, 'last');
            for i = 1:nb
                B = m.Bodies(i);
                if isempty(B.F)
                    set(h(i), 'XData', x.r(1,i,k), 'YData', x.r(2,i,k), 'ZData', x.r(3,i,k));
                else
                    set(h(i), 'Vertices', B.Xl*x.R(:,:,i,k)' + x.r(:,i,k)');
                end
            end
            title(sprintf('%s  t = %.3f s', opts.Title, t(k)));
            drawnow limitrate
        end
    end
end

methods (Access = private)
    function addJoint(m, name, iA, iB, pA, RA, pB, RB, free, k, c)
    % Joint between reference points iA (0 = ground) and iB: frame (pA, RA)
    % on A and (pB, RB) on B, in world coordinates at the initial
    % configuration (body frames are world-aligned there).
        bA = 0;  a = pA;  QA = RA;
        if iA > 0, bA = m.Refs(iA).body;  a = pA - m.Bodies(bA).c0; end
        bB = m.Refs(iB).body;
        assert(bA ~= bB, 'fembd:Model', '%s connects a body to itself.', name);
        kk = zeros(6, 1);  cc = zeros(6, 1);
        kk(free) = k(:).*ones(nnz(free), 1);  cc(free) = c(:).*ones(nnz(free), 1);
        m.Joints(end+1) = struct('name', name, 'A', bA, 'B', bB, 'a', a, ...
            'b', pB - m.Bodies(bB).c0, 'QA', QA, 'QB', RB, 'free', logical(free), 'k', kk, 'c', cc);
    end

    function [s, Js, g, Gr] = jointKin(m, j, q)
    % Joint coordinates s = [d; phi] (translation in the A frame, rotation
    % vector of B relative to A) with Jacobian Js, and the rotation
    % constraints g (dot products of frame axes) with Jacobian Gr.
        J = m.Joints(j);  n = 6*numel(m.Bodies);
        [rA, RA] = deal(zeros(3, 1), eye(3));
        if J.A > 0, rA = q.r(:, J.A);  RA = q.R(:,:,J.A); end
        rB = q.r(:, J.B);  RB = q.R(:,:,J.B);
        D = rB + RB*J.b - rA - RA*J.a;
        Rrel = J.QA'*RA'*RB*J.QB;
        phi = logso3(Rrel);  Ti = dexpInv(phi);
        s = [J.QA'*RA'*D; phi];
        Js = zeros(6, n);
        cB = 6*J.B + (-5:0);
        Js(1:3, cB(1:3)) = J.QA'*RA';
        Js(1:3, cB(4:6)) = -J.QA'*RA'*RB*hat(J.b);
        Js(4:6, cB(4:6)) = Ti*J.QB';
        if J.A > 0
            cA = 6*J.A + (-5:0);
            Js(1:3, cA(1:3)) = -J.QA'*RA';
            Js(1:3, cA(4:6)) = J.QA'*(hat(RA'*D) + hat(J.a));
            Js(4:6, cA(4:6)) = -Ti*Rrel'*J.QA';
        end
        fr = J.free(4:6);
        switch nnz(fr)
            case 0, pairs = [2 3; 3 1; 1 2];
            case 1, k = find(fr);  i = setdiff(1:3, k);  pairs = [i(:), [k; k]];
            case 2, k = find(~fr);  pairs = circshift([1 2 3], -k);  pairs = pairs(1:2);
            otherwise, pairs = zeros(0, 2);
        end
        g = zeros(size(pairs, 1), 1);  Gr = zeros(size(pairs, 1), n);
        for p = 1:size(pairs, 1)
            ua = J.QA(:, pairs(p,1));  wb = J.QB(:, pairs(p,2));
            u = RA*ua;  w = RB*wb;
            g(p) = u'*w;
            Gr(p, cB(4:6)) = -cross(RB'*u, wb)';
            if J.A > 0, Gr(p, cA(4:6)) = -cross(RA'*w, ua)'; end
        end
    end

    function [v, Jv] = coord(m, kind, idx, q)
    % Coordinates of a joint (its free coordinates) or reference point
    % ([displacement; rotation vector], absolute) and their Jacobian.
        if strcmp(kind, 'joint')
            [s, Js] = m.jointKin(idx, q);  f = m.Assembly.Joints(idx).free;
            v = s(f);  Jv = Js(f, :);
        else
            R = m.Refs(idx);  i = R.body;  Ri = q.R(:,:,i);  c = 6*i + (-5:0);
            phi = logso3(Ri);
            v = [q.r(:, i) + Ri*R.s - R.x0; phi];
            Jv = zeros(6, 6*numel(m.Bodies));
            Jv(1:3, c(1:3)) = eye(3);  Jv(1:3, c(4:6)) = -Ri*hat(R.s);
            Jv(4:6, c(4:6)) = dexpInv(phi);
        end
    end

    function [d, Jd] = stretch(m, T, q)
    % Transmission stretch delta = sum c_i s_i (minus its initial value).
        d = -T.d0;  Jd = zeros(1, 6*numel(m.Bodies));
        for i = 1:numel(T.terms)
            [v, Jv] = m.coord(T.terms(i).kind, T.terms(i).idx, q);
            d = d + T.coef(i)*v(T.terms(i).dir);  Jd = Jd + T.coef(i)*Jv(T.terms(i).dir, :);
        end
    end

    function [g, G] = constraints(m, q, reduce)
    % Constraint values g(q) = 0 and their Jacobian G (rows m.keep).
        if nargin < 3, reduce = true; end
        n = 6*numel(m.Bodies);  g = zeros(0, 1);  G = zeros(0, n);
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

    function [F, Kt, Ct] = forces(m, q, nu, u)
    % Generalised forces F (gravity, gyroscopic, springs, dampers, inputs)
    % and approximate tangents Kt = -dF/dq, Ct = -dF/dnu.
        nb = numel(m.Bodies);  n = 6*nb;
        F = zeros(n, 1);  Kt = zeros(n);  Ct = zeros(n);
        for i = 1:nb
            c = 6*i + (-5:0);  w = nu(c(4:6));  J = m.Bodies(i).J;
            F(c(1:3)) = m.Bodies(i).m*m.Gravity(:);
            F(c(4:6)) = -cross(w, J*w);
            Ct(c(4:6), c(4:6)) = hat(w)*J - hat(J*w);
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
            if strcmp(s.kind, 'joint')
                [~, Jv] = m.coord('joint', s.idx, q);  F = F + Jv(s.dir, :)'*u(i);
            else
                R = m.Refs(s.idx);  b = R.body;  c = 6*b + (-5:0);  Rb = q.R(:,:,b);
                e = zeros(3, 1);  e(mod(s.dir - 1, 3) + 1) = u(i);
                if s.dir <= 3                    % force at the reference point
                    F(c(1:3)) = F(c(1:3)) + e;
                    F(c(4:6)) = F(c(4:6)) + cross(R.s, Rb'*e);
                else                             % moment (world axes)
                    F(c(4:6)) = F(c(4:6)) + Rb'*e;
                end
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
    % Input/output selections -> signals (kind, index, direction).
        a = m.Assembly;
        specs = struct('kind', {}, 'idx', {}, 'dir', {}, 'type', {});  names = {};
        for s = list
            assert(isempty(s.nodes), 'fembd:Model', ...
                '%s: node-set signals are not supported by the nonlinear model yet.', s.Name);
            [kind, k] = a.lookup(s.id, true);
            dirs = s.dirs;
            if isempty(dirs)
                dirs = 1:6;  if strcmp(kind, 'joint'), dirs = 1:nnz(a.Joints(k).free); end
            end
            type = 'displacement';  if isfield(s, 'type'), type = s.type; end
            for d = dirs(:)'
                specs(end+1) = struct('kind', kind, 'idx', k, 'dir', d, 'type', type); %#ok<AGROW>
                names{end+1} = s.Name;                                              %#ok<AGROW>
                if numel(dirs) > 1, names{end} = sprintf('%s_%d', s.Name, d); end
            end
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
end

methods (Static, Access = private)
    function q = retract(q, dq)
    % q o exp(dq): positions add (world), rotations update on the right (body).
        for i = 1:size(q.r, 2)
            c = 6*i + (-5:0);
            q.r(:, i) = q.r(:, i) + dq(c(1:3));
            q.R(:,:,i) = q.R(:,:,i)*expso3(dq(c(4:6)));
        end
    end
end
end

function [m, c, J] = massProps(b)
% Mass, centre of mass and inertia (about c) of a rigid body.
if ~isempty(b.MassProps)
    p = b.MassProps;  m = p.m;  c = p.c(:)';  J = p.J;  return
end
[~, M] = b.matrices();  Ms = M(1:3:end, 1:3:end);
m = full(sum(Ms(:)));
c = full(sum(Ms, 1)*b.X)/m;
Xc = b.X - c;  S = full(Xc'*Ms*Xc);
J = trace(S)*eye(3) - S;  J = (J + J')/2;
end

function S = hat(v)
S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function R = expso3(p)
% Rotation matrix exp(hat(p)) (Rodrigues).
th = norm(p);  P = hat(p);
if th < 1e-8, R = eye(3) + P + P*P/2; return, end
R = eye(3) + sin(th)/th*P + (1 - cos(th))/th^2*(P*P);
end

function p = logso3(R)
% Rotation vector of R (angle in [0, pi]).
s = [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)]/2;
c = (trace(R) - 1)/2;  sn = norm(s);  th = atan2(sn, c);
if th < 1e-8
    p = s;
elseif pi - th > 1e-6
    p = th/sn*s;
else                                         % near pi: axis from R + I
    B = (R + eye(3))/2;  [~, i] = max(diag(B));
    n = B(:, i)/sqrt(B(i, i));
    if n'*s < 0, n = -n; end
    p = th*n;
end
end

function Ti = dexpInv(p)
% Inverse of the right Jacobian of SO(3): dp/dt = Ti*omega for R = exp(hat(p)),
% dR/dt = R*hat(omega).
th = norm(p);  P = hat(p);
if th < 1e-6
    Ti = eye(3) + P/2 + (P*P)/12;
else
    Ti = eye(3) + P/2 + (1/th^2 - (1 + cos(th))/(2*th*sin(th)))*(P*P);
end
end
