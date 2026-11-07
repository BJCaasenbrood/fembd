classdef Assembly < handle
%ASSEMBLY Elastic and rigid bodies with ties, couplings, joints and supports.
%   a = fem.Assembly();  a.add(body);  a.tie(master, slave);  a.fix(...);
%   r = a.couple(body, sel);  a.joint(A, B, 'revolute', Axis = ...);
%   [K, M, T] = a.system();   % reduced matrices, full displacement u = T*q
%
%   Bodies are referenced by index, name or handle; reference points (rigid
%   bodies, couplings) and joints by name. Global dofs are, in order:
%     node dofs of all bodies [ux1 uy1 uz1 ux2 ...], in the order added
%     6 dofs [ux uy uz rx ry rz] per reference point, in the order created
%     the released dofs of each joint, in the order created
%   Nodes of a rigid body or coupling follow their reference point:
%   u = u_ref + r x (x - x_ref).

properties
    Bodies = {}
    Links  = struct('slave', {}, 'master', {}, 'N', {}, 'n', {}, 'A', {}, 'k', {})
    Refs   = struct('Name', {}, 'body', {}, 'nodes', {}, 'x', {})   % reference points
    Joints = struct('Name', {}, 'A', {}, 'B', {}, 'p', {}, 'R', {}, 'free', {}, 'k', {}, 'value', {})
    Fixes  = struct('id', {}, 'nodes', {}, 'dirs', {})              % supports (u = 0)
end

properties (Dependent)
    nDof
    X       % all nodes, stacked
end

methods
    function b = add(a, b)
    %ADD Append a fem.Body. A rigid body also gets a reference point with
    %   the body's name.
        a.Bodies{end+1} = b;
        if b.Rigid
            a.Refs(end+1) = struct('Name', b.Name, 'body', b, 'nodes', [], 'x', []);
        end
    end

    function n = get.nDof(a)
        off = a.offsets;
        n = 3*off(end) + 6*numel(a.Refs) + sum(arrayfun(@(J) nnz(J.free), a.Joints));
    end
    function X = get.X(a),    X = cell2mat(cellfun(@(b) b.X, a.Bodies(:), 'UniformOutput', false)); end

    function d = dof(a, id, sel, dirs)
    %DOF Global dofs.
    %   a.dof(body, sel, dirs)  nodes selected by sel(X) (default all), dirs
    %                           subset of 1:3 (default 1:3)
    %   a.dof(ref, [], dirs)    reference point (rigid body or coupling),
    %                           dirs subset of 1:6 [ux uy uz rx ry rz]
    %   a.dof(joint, [], dirs)  released joint coordinates
        if nargin < 3, sel = []; end
        if nargin < 4, dirs = []; end
        [kind, k] = a.lookup(id, isempty(sel));
        switch kind
            case 'joint'
                d = a.jointDofs(k);
            case 'ref'
                d = a.refDofs(k);
            otherwise
                b = a.Bodies{k};  off = a.offsets;
                if isempty(sel), sel = (1:b.nNode)'; end
                if isempty(dirs), dirs = 1:3; end
                d = reshape(3*(off(k) + b.nodes(sel))' - 3 + dirs(:), [], 1);
                return
        end
        if ~isempty(dirs), d = d(dirs(:)); end
    end

    function fix(a, id, sel, dirs)
    %FIX Support (u = 0) the dofs selected as in a.dof: nodes of an elastic
    %   body, or a.fix(ref, [], dirs) / a.fix(joint, [], dirs).
        if nargin < 3, sel = []; end
        if nargin < 4, dirs = []; end
        [kind, k] = a.lookup(id, isempty(sel));
        nodes = [];
        if strcmp(kind, 'body')
            b = a.Bodies{k};  id = b;
            assert(~b.Rigid, 'fem:fix', ...
                'Rigid bodies are supported via their reference dofs: a.fix(body, [], dirs).');
            nodes = (1:b.nNode)';
            if ~isempty(sel), nodes = b.nodes(sel); end
            if isempty(nodes), warning('fem:fix', 'No nodes selected.'); return, end
        end
        a.Fixes(end+1) = struct('id', {id}, 'nodes', nodes, 'dirs', dirs);
    end

    function name = couple(a, body, sel, opts)
    %COUPLE Rigid coupling (RBE2): the nodes of an elastic body selected by
    %   sel follow a new 6-dof reference point, e.g. a bore for a joint.
    %   Returns the reference name. Options: Name, Point (default: centroid
    %   of the selected nodes).
        arguments
            a
            body
            sel
            opts.Name (1,:) char = ''
            opts.Point double = []
        end
        b = a.body(body);
        assert(~b.Rigid, 'fem:couple', 'Rigid bodies already have a reference point.');
        nodes = b.nodes(sel);
        assert(~isempty(nodes), 'fem:couple', 'No nodes selected.');
        x = opts.Point;
        if isempty(x), x = mean(b.X(nodes,:), 1); end
        name = opts.Name;
        if isempty(name), name = sprintf('%s_ref%d', b.Name, numel(a.Refs) + 1); end
        a.Refs(end+1) = struct('Name', name, 'body', b, 'nodes', nodes, 'x', x(:)');
    end

    function name = joint(a, A, B, type, opts)
    %JOINT Kinematic joint between reference points A and B (rigid bodies,
    %   couplings, or 'ground' for A). B follows A, except for the released
    %   dofs, which become joint coordinates (a.dof(name)). Joints that close
    %   a kinematic loop are handled as cut joints (see constraints).
    %
    %   type      'fixed' (default), 'revolute', 'prismatic', 'cylindrical',
    %             'spherical', 'planar', 'free' (bushing, use Stiffness)
    %   Point     joint location (default: reference point of B)
    %   Axis      joint axis = local z (default [0 0 1]); revolute and
    %             prismatic act about/along it, planar moves normal to it
    %   Frame     3x3 [ex ey ez] instead of Axis
    %   Free      1x6 logical mask [ux uy uz rx ry rz] in the joint frame,
    %             overrides type
    %   Stiffness stiffness of the released dofs (scalar or per dof)
    %   Name      default 'joint<k>'
        arguments
            a
            A
            B
            type (1,:) char {mustBeMember(type, {'fixed', 'revolute', 'prismatic', ...
                'cylindrical', 'spherical', 'planar', 'free'})} = 'fixed'
            opts.Point double = []
            opts.Axis (1,3) double = [0 0 1]
            opts.Frame double = []
            opts.Free = []
            opts.Stiffness double = 0
            opts.Name (1,:) char = ''
        end
        iA = a.refIndex(A);  iB = a.refIndex(B);
        assert(iB > 0 && iA ~= iB, 'fem:joint', 'B must be a reference point other than A.');
        mask = struct('fixed', [0 0 0 0 0 0], 'revolute', [0 0 0 0 0 1], ...
            'prismatic', [0 0 1 0 0 0], 'cylindrical', [0 0 1 0 0 1], ...
            'spherical', [0 0 0 1 1 1], 'planar', [1 1 0 0 0 1], 'free', [1 1 1 1 1 1]);
        free = logical(mask.(type));
        if ~isempty(opts.Free), free = logical(opts.Free(:)'); end
        R = opts.Frame;
        if isempty(R)
            ez = opts.Axis/norm(opts.Axis);
            h = [1 0 0];  if abs(ez(1)) > 0.9, h = [0 1 0]; end
            ex = h - dot(h, ez)*ez;  ex = ex/norm(ex);
            R = [ex; cross(ez, ex); ez]';
        end
        p = opts.Point;
        if isempty(p), p = a.refPoint(iB); end
        name = opts.Name;
        if isempty(name), name = sprintf('joint%d', numel(a.Joints) + 1); end
        a.Joints(end+1) = struct('Name', name, 'A', iA, 'B', iB, 'p', p(:)', 'R', R, ...
            'free', free, 'k', opts.Stiffness(:).*ones(nnz(free), 1), 'value', zeros(nnz(free), 1));
    end

    function ns = tie(a, master, slave, opts)
    %TIE Connect the surface nodes of slave to the surface of master.
    %   Slave surface nodes within Tol of the master surface are coupled to
    %   the closest master point, interpolated with the master face shape
    %   functions. Without Stiffness the coupling is a rigid (bonded) tie,
    %   eliminated in T. With Stiffness = [kn kt] (per unit area) it is an
    %   elastic interface added to K. Use the finer mesh as slave. In a rigid
    %   tie with a rigid body, the rigid body is always made the master.
    %
    %   Options: Tol (default 0.1 x median face edge of the finer surface), Stiffness,
    %            Select (selector @(X) on slave nodes to restrict the tie)
        arguments
            a
            master
            slave
            opts.Tol double = []
            opts.Stiffness double = []
            opts.Select = []
        end
        [bm, offm] = a.body(master);  [bs, offs] = a.body(slave);
        assert(~isempty(bm.F) && ~isempty(bs.F), 'fem:tie', ...
            'A point mass has no surface; connect it with a.joint.');
        if bs.Rigid && isempty(opts.Stiffness)
            assert(~bm.Rigid, 'fem:tie', ['Two rigid bodies cannot be tied rigidly; ' ...
                'make them one body or use Stiffness.']);
            [bm, offm, bs, offs] = deal(bs, offs, bm, offm);
        end
        Fm = bm.F;
        A = bm.X(Fm(:,1),:);  B = bm.X(Fm(:,2),:);  C = bm.X(Fm(:,3),:);
        tol = opts.Tol;
        if isempty(tol)
            edge = @(b, F) median(vecnorm(b.X(F(:,[2 3 1]),:) - b.X(F(:,1:3),:), 2, 2));
            tol = 0.1*min(edge(bm, Fm), edge(bs, bs.F));
        end

        s = unique(bs.F(:));
        if ~isempty(opts.Select), s = s(opts.Select(bs.X(s,:))); end
        lo = min(min(A, B), C) - tol;  hi = max(max(A, B), C) + tol;
        s = s(all(bs.X(s,:) >= min(lo) & bs.X(s,:) <= max(hi), 2));

        d = inf(numel(s), 1);  f = zeros(numel(s), 1);  l = zeros(numel(s), 3);
        for i = 1:numel(s)
            p = bs.X(s(i),:);
            c = find(all(p >= lo & p <= hi, 2));
            if isempty(c), continue, end
            [dc, lc] = closestOnTri(p, A(c,:), B(c,:), C(c,:));
            [d(i), k] = min(dc);  f(i) = c(k);  l(i,:) = lc(k,:);
        end
        on = d <= tol;  s = s(on);  f = f(on);  l = l(on,:);
        ns = numel(s);
        if ns == 0, warning('fem:tie', 'No slave nodes within tolerance %g.', tol); return, end

        n = cross(B(f,:) - A(f,:), C(f,:) - A(f,:), 2);
        n = n./vecnorm(n, 2, 2);

        % nodal areas of the slave interface (row-sum / HRZ face weights) from
        % slave faces that touch a tied node and face the master surface
        Fs = bs.F;
        nf = cross(bs.X(Fs(:,2),:) - bs.X(Fs(:,1),:), bs.X(Fs(:,3),:) - bs.X(Fs(:,1),:), 2);
        nm = zeros(bs.nNode, 3);  nm(s,:) = n;
        nm = reshape(sum(reshape(nm(Fs,:), [size(Fs) 3]), 2), [], 3);
        on = dot(nf, nm, 2) < -0.5*vecnorm(nf, 2, 2).*vecnorm(nm, 2, 2);
        Fs = Fs(on,:);  ar = vecnorm(nf(on,:), 2, 2)/2;
        w  = [1 1 1]/3;
        if size(Fs, 2) == 6, w = [3 3 3 16 16 16]/57; end
        nodeA = accumarray(Fs(:), reshape(ar.*w, [], 1), [bs.nNode 1]);

        a.Links(end+1) = struct('slave', offs + s, 'master', offm + Fm(f,:), ...
            'N', triShape(l, size(Fm, 2)), 'n', n, 'A', nodeA(s), 'k', opts.Stiffness);
    end

    function [K, M] = assemble(a, opts)
    %ASSEMBLE Unconstrained global K and M, including elastic interfaces.
    %   Option Lumped = true gives an HRZ-lumped (diagonal) mass matrix.
        arguments
            a
            opts.Lumped (1,1) logical = false
        end
        [K, M] = cellfun(@(b) b.matrices(Lumped = opts.Lumped), a.Bodies, 'UniformOutput', false);
        n = a.nDof;  off = a.offsets;  nx = n - 3*off(end);
        K = blkdiag(K{:}, sparse(nx, nx));  M = blkdiag(M{:}, sparse(nx, nx));
        for k = 1:numel(a.Refs)                           % mesh-free rigid bodies
            b = a.Refs(k).body;
            if isempty(a.Refs(k).nodes) && ~isempty(b.MassProps)
                d = a.refDofs(k);
                M(d, d) = M(d, d) + b.rigidMass();
            end
        end
        for k = 1:numel(a.Joints)                         % joint springs
            d = a.jointDofs(k);
            K = K + sparse(d, d, a.Joints(k).k, n, n);
        end
        for L = a.Links(~arrayfun(@(L) isempty(L.k), a.Links))
            ns = numel(L.slave);  kn = L.k(1);  kt = L.k(end);
            nn = permute(L.n, [2 3 1]).*permute(L.n, [3 2 1]);
            Db = reshape(L.A, 1, 1, []).*(kt*eye(3) + (kn - kt)*nn);
            [p, q, i] = ndgrid(1:3, 1:3, 1:ns);
            D = sparse(3*(i(:) - 1) + p(:), 3*(i(:) - 1) + q(:), Db(:));
            G = a.gap(L);
            K = K + G'*D*G;
        end
    end

    function [T, I] = constraints(a)
    %CONSTRAINTS Map u = T*q from independent dofs q to all dofs u, which
    %   eliminates rigid bodies, couplings, joints, rigid ties (slave dofs)
    %   and supports. I lists the independent dofs: u(I) = q.
    %   Joints in a tree are eliminated master-slave (B follows A). A joint
    %   that closes a loop (or whose B side already follows another joint or
    %   is supported) is a cut joint: its constrained directions give rows
    %   C*q = 0, removed by rank-revealing QR coordinate partitioning, which
    %   also drops redundant rows (e.g. planar loops in 3D).
        n = a.nDof;  off = a.offsets;
        S = zeros(0, 1);  R = sparse(0, n);
        for k = 1:numel(a.Refs)                  % nodes following a reference point
            [g, X] = a.refNodes(k);  r = a.refDofs(k);
            nn = numel(g);  d = X - a.refPoint(k);  j = 3*(0:nn-1)';
            S = [S; reshape(3*g' - [2; 1; 0], [], 1)];                    %#ok<AGROW>
            R = [R; sparse([j+1; j+2; j+3; j+1; j+1; j+2; j+2; j+3; j+3], ...
                           repelem(r([1 2 3 5 6 4 6 4 5]), nn), ...
                           [ones(3*nn, 1); d(:,3); -d(:,2); -d(:,3); d(:,1); d(:,2); -d(:,1)], ...
                           3*nn, n)];                                     %#ok<AGROW>
        end
        skew = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        fixed = a.fixedDofs();  cut = a.cutJoints(fixed);  Cc = sparse(0, n);
        for k = 1:numel(a.Joints)
            J = a.Joints(k);  E = eye(6);  Rh = blkdiag(J.R, J.R);
            Pinv = [eye(3), skew(J.p - a.refPoint(J.B)); zeros(3), eye(3)];
            G = sparse(6, n);                    % relative motion at p: P_B q_B - P_A q_A
            G(:, a.refDofs(J.B)) = inv(Pinv);
            if J.A > 0
                G(:, a.refDofs(J.A)) = -[eye(3), -skew(J.p - a.refPoint(J.A)); zeros(3), eye(3)];
            end
            if ~cut(k)                           % q_B = P_B^-1 (P_A q_A + R E s)
                Rk = -Pinv*G;  Rk(:, a.refDofs(J.B)) = 0;
                Rk(:, a.jointDofs(k)) = Pinv*Rh*E(:, J.free);
                S = [S; a.refDofs(J.B)];  R = [R; Rk];                    %#ok<AGROW>
            else                                 % s = (R E)' G q,  (R E_c)' G q = 0
                S = [S; a.jointDofs(k)];  R = [R; (Rh*E(:, J.free))'*G];  %#ok<AGROW>
                Cc = [Cc; (Rh*E(:, ~J.free))'*G];                         %#ok<AGROW>
            end
        end
        for L = a.Links(arrayfun(@(L) isempty(L.k), a.Links))
            s = reshape(3*L.slave' - [2; 1; 0], [], 1);
            S = [S; s];                                              %#ok<AGROW>
            R = [R; sparse(1:numel(s), s, 1, numel(s), n) - a.gap(L)]; %#ok<AGROW>
        end
        [S, i] = unique(S, 'stable');  R = R(i,:);
        keep = ~ismember(S, fixed);  S = S(keep);  R = R(keep,:);

        % resolve chains (a slave that is master elsewhere)
        mask = sparse(S, S, 1, n, n);
        for it = 1:100
            C = R(:, S);
            if nnz(C) == 0, break, end
            R = R - R*mask + C*R;
        end
        assert(nnz(R(:, S)) == 0, 'fem:constraints', 'Cyclic constraints (e.g. a closed loop of ties or joints).');

        I = setdiff((1:n)', [S; fixed]);
        [i, j, v] = find(R(:, I));
        T = sparse(I, 1:numel(I), 1, n, numel(I)) + sparse(S(i), j, v, n, numel(I));

        if ~isempty(Cc)                          % closed loops: C q = 0
            pr = inf(numel(I), 1);               % eliminate non-joint dofs first,
            for k = find(~cut)                   % then later joints before earlier
                pr(ismember(I, a.jointDofs(k))) = k;
            end
            [~, order] = sort(pr, 'descend');
            [T2, keep] = partition(Cc*T, order);
            T = T*T2;  I = I(keep);
        end
    end

    function [K, M, T] = system(a, opts)
    %SYSTEM Reduced matrices K = T'*K*T, M = T'*M*T and the map u = T*q.
    %   Option Lumped = true lumps the body mass matrices (HRZ). M stays
    %   diagonal unless rigid ties redistribute slave masses to masters.
        arguments
            a
            opts.Lumped (1,1) logical = false
        end
        [K, M] = a.assemble(Lumped = opts.Lumped);
        T = a.constraints();
        K = T'*K*T;  K = (K + K')/2;
        M = T'*M*T;  M = (M + M')/2;
    end

    function pose(a, name, value)
    %POSE Move a joint to a finite position and rebuild about it.
    %   a.pose(joint, value) sets the joint coordinates (rad for rotations,
    %   m for translations, in the joint frame) relative to the configuration
    %   in which the joint was created; pose(joint, 0) restores it. All parts
    %   on the B side (bodies, couplings, downstream joints and tied bodies)
    %   move rigidly. The bodies themselves are moved (they are handles).
    %   K and M from a.system() are then linearised about the new pose. For
    %   several rotations (spherical joint) value is a rotation vector.
        [kind, k] = a.lookup(name, true);
        assert(strcmp(kind, 'joint'), 'fem:pose', '%s is not a joint.', name);
        J = a.Joints(k);  value = value(:);
        assert(numel(value) == nnz(J.free), 'fem:pose', ...
            '%s has %d joint coordinate(s).', J.Name, nnz(J.free));
        [Rc, tc] = jointMotion(J, J.value);  [Rn, tn] = jointMotion(J, value);
        Rd = Rn*Rc';  td = J.p(:) + tn - Rd*(J.p(:) + tc);  % x -> Rd*x + td
        moving = a.subtree(k);
        for i = find(moving), a.Bodies{i}.transform(Rd, td); end
        for r = 1:numel(a.Refs)
            if ~isempty(a.Refs(r).x) && moving(a.refBody(r))
                a.Refs(r).x = (Rd*a.Refs(r).x(:) + td)';
            end
        end
        for j = 1:numel(a.Joints)
            Jj = a.Joints(j);
            if j ~= k && Jj.A > 0 && moving(a.refBody(Jj.A))
                a.Joints(j).p = (Rd*Jj.p(:) + td)';  a.Joints(j).R = Rd*Jj.R;
            end
        end
        for l = 1:numel(a.Links)
            if moving(a.nodeBody(a.Links(l).master(1)))
                a.Links(l).n = a.Links(l).n*Rd';
            end
        end
        a.Joints(k).value = value;
    end

    function j = qdof(a, d)
    %QDOF Indices into q of global dofs d (e.g. from a.dof), so u(d) = q(j).
    %   The dofs must be independent: not supported, not dependent (tie
    %   slaves, nodes of rigid bodies/couplings, B sides of joints), and not
    %   eliminated by a closed loop.
        [~, I] = a.constraints();
        [ok, j] = ismember(d, I);
        assert(all(ok), 'fem:qdof', ['Some dofs are supported, dependent, or ' ...
            'eliminated by a closed loop.']);
    end

    function plot(a, u, scale)
    %PLOT Draw the boundary of all bodies, optionally deformed by u (nDof x 1).
        if nargin < 2 || isempty(u), u = zeros(a.nDof, 1); end
        if nargin < 3, scale = 1; end
        off = a.offsets;  U = reshape(u(1:3*off(end)), 3, [])';  c = lines(numel(a.Bodies));
        held = ishold;  hold on
        for i = 1:numel(a.Bodies)
            b = a.Bodies{i};  Ui = U(off(i)+1:off(i+1), :);  F = b.F;
            if isempty(F)                                 % point mass: marker
                [~, k] = a.lookup(b, true);  x = b.Ref;  d = a.refDofs(k);
                if any(u), x = x + scale*u(d(1:3))'; end
                plot3(x(1), x(2), x(3), 'ko', 'MarkerFaceColor', c(i,:), 'MarkerSize', 8);
                continue
            end
            if size(F, 2) == 6
                F = [F(:,[1 4 6]); F(:,[4 2 5]); F(:,[6 5 3]); F(:,[4 5 6])];
            end
            if any(u)
                patch('Faces', F, 'Vertices', b.X + scale*Ui, 'FaceVertexCData', ...
                      vecnorm(Ui, 2, 2), 'FaceColor', 'interp', 'EdgeAlpha', 0.2);
            else
                patch('Faces', F, 'Vertices', b.X, 'FaceColor', c(i,:), 'EdgeAlpha', 0.2);
            end
        end
        if ~held, hold off, end
        axis equal;  view(3);
    end
end

methods (Access = private)
    function off = offsets(a)
        off = cumsum([0, cellfun(@(b) b.nNode, a.Bodies)]);
    end

    function [b, off, i] = body(a, id)
        if isa(id, 'fem.Body')
            i = find(cellfun(@(b) b == id, a.Bodies));
        elseif ischar(id) || isstring(id)
            i = find(strcmp(cellfun(@(b) b.Name, a.Bodies, 'UniformOutput', false), id));
        else
            i = id;
        end
        assert(isscalar(i), 'fem:Assembly', 'Body not found or ambiguous.');
        b = a.Bodies{i};  off = a.offsets;  off = off(i);
    end

    function cut = cutJoints(a, fixed)
    % Joints that cannot be eliminated master-slave: their ends are already
    % connected through earlier joints (loop), their B side already follows
    % another joint, or their B side is supported. Ground is node 0.
        nr = numel(a.Refs);  root = 0:nr;  isB = false(1, nr);
        cut = false(1, numel(a.Joints));
        for k = 1:numel(a.Joints)
            J = a.Joints(k);  rA = J.A;  rB = J.B;
            while root(rA+1) ~= rA, rA = root(rA+1); end
            while root(rB+1) ~= rB, rB = root(rB+1); end
            cut(k) = isB(J.B) || rA == rB || any(ismember(a.refDofs(J.B), fixed));
            if ~cut(k), isB(J.B) = true;  root(rB+1) = rA; end
        end
    end

    function moving = subtree(a, k)
    % Bodies on the B side of joint k: connected through ties or child joints.
        moving = false(1, numel(a.Bodies));  J = a.Joints(k);
        moving(a.refBody(J.B)) = true;
        E = zeros(numel(a.Links), 2);
        for l = 1:numel(a.Links)
            E(l,:) = [a.nodeBody(a.Links(l).master(1)), a.nodeBody(a.Links(l).slave(1))];
        end
        others = setdiff(1:numel(a.Joints), k);
        grown = true;
        while grown
            old = moving;
            moving(E(any(moving(E), 2), :)) = true;
            for j = others
                Jj = a.Joints(j);
                if Jj.A > 0 && moving(a.refBody(Jj.A)), moving(a.refBody(Jj.B)) = true; end
            end
            grown = any(moving ~= old);
        end
        loop = (J.A > 0 && moving(a.refBody(J.A))) || any(arrayfun(@(j) ...
            moving(a.refBody(a.Joints(j).B)) && (a.Joints(j).A == 0 || ...
            ~moving(a.refBody(a.Joints(j).A))), others));
        for f = a.Fixes
            [kind, i] = a.lookup(f.id, isempty(f.nodes));
            if strcmp(kind, 'ref'), i = a.refBody(i); end
            loop = loop || (~strcmp(kind, 'joint') && moving(i));
        end
        assert(~loop, 'fem:pose', ['%s: the moving side is also attached elsewhere ' ...
            '(closed loop or support), so it cannot be posed.'], J.Name);
    end

    function i = refBody(a, r)
        i = find(cellfun(@(b) b == a.Refs(r).body, a.Bodies));
    end

    function i = nodeBody(a, g)
        off = a.offsets;  i = find(g > off(1:end-1), 1, 'last');
    end

    function [kind, k] = lookup(a, id, whole)
    % Resolve id to a 'joint', 'ref' or 'body' index. Joints and reference
    % points take precedence when no node selection is given (whole).
        if whole && (ischar(id) || isstring(id))
            k = find(strcmp({a.Joints.Name}, id));  kind = 'joint';
            if isscalar(k), return, end
            k = find(strcmp({a.Refs.Name}, id));    kind = 'ref';
            if isscalar(k), return, end
        end
        [b, ~, k] = a.body(id);  kind = 'body';
        if whole && b.Rigid
            k = find(arrayfun(@(r) r.body == b && isempty(r.nodes), a.Refs));  kind = 'ref';
        end
    end

    function k = refIndex(a, id)
    % Reference point index of id; 0 for ground.
        k = 0;
        if isempty(id) || (ischar(id) || isstring(id)) && strcmpi(id, 'ground'), return, end
        [kind, k] = a.lookup(id, true);
        assert(strcmp(kind, 'ref'), 'fem:joint', ...
            'Joint ends must be rigid bodies, couplings (a.couple) or ''ground''.');
    end

    function [g, X] = refNodes(a, k)
    % Global indices and coordinates of the nodes following reference k.
        r = a.Refs(k);  i = find(cellfun(@(b) b == r.body, a.Bodies));  off = a.offsets;
        nodes = r.nodes;
        if isempty(nodes), nodes = (1:r.body.nNode)'; end
        g = off(i) + nodes;  X = r.body.X(nodes,:);
    end

    function x = refPoint(a, k)
        x = a.Refs(k).x;
        if isempty(x), x = a.Refs(k).body.Ref; end
    end

    function r = refDofs(a, k)
        off = a.offsets;
        r = 3*off(end) + 6*(k - 1) + (1:6)';
    end

    function d = jointDofs(a, k)
        nf = arrayfun(@(J) nnz(J.free), a.Joints);  off = a.offsets;
        d = 3*off(end) + 6*numel(a.Refs) + sum(nf(1:k-1)) + (1:nf(k))';
    end

    function d = fixedDofs(a)
    % Global dofs of all supports.
        d = arrayfun(@(f) a.dof(f.id, f.nodes, f.dirs), a.Fixes, 'UniformOutput', false);
        d = unique(vertcat(zeros(0, 1), d{:}));
    end

    function G = gap(a, L)
    % Sparse gap operator g = G*u = u_slave - sum_j N_j u_master,j.
        d = reshape(1:3, 1, 1, 3);
        nodes = [L.slave, L.master];
        rows = 3*(0:numel(L.slave) - 1)' + d + 0*nodes;
        cols = 3*nodes - 3 + d;
        vals = [ones(numel(L.slave), 1), -L.N] + 0*d;
        G = sparse(rows(:), cols(:), vals(:), 3*numel(L.slave), a.nDof);
    end
end
end

function [R, t] = jointMotion(J, v)
% Finite rigid motion of a joint about its point: x -> R*(x - p) + p + t.
s = zeros(6, 1);  s(J.free) = v;
t = J.R*s(1:3);  w = J.R*s(4:6);
R = expm([0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0]);
end

function [T2, keep] = partition(C, order)
% Coordinate partitioning of C*q = 0: q = T2*q(keep). Dependent columns are
% picked greedily in the given priority order, each only if it raises the
% rank, so redundant rows (e.g. planar loops in 3D) are dropped. The
% dependent dofs then follow by least squares (exact for consistent C).
m = size(C, 2);
C = C./max(vecnorm(C, 2, 2), eps);              % unit rows
cols = order(any(C(:, order), 1));  cols = cols(:)';
Cn = full(C(:, cols));  Cn = Cn./vecnorm(Cn);   % unit columns for rank tests
tol = 1e-9;  target = rank(Cn, tol);  sel = [];
for c = 1:numel(cols)
    if numel(sel) == target, break, end
    if rank(Cn(:, [sel c]), tol) > numel(sel), sel(end+1) = c; end %#ok<AGROW>
end
dep = cols(sel);  ind = cols(setdiff(1:numel(cols), sel));
X = -full(C(:, dep))\full(C(:, ind));
keep = setdiff(1:m, dep);
[~, pos] = ismember(ind, keep);
T2 = sparse(keep, 1:numel(keep), 1, m, numel(keep));
[i, j, v] = find(X);
T2 = T2 + sparse(dep(i), pos(j), v, m, numel(keep));
end
