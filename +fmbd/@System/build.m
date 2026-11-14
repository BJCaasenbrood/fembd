function build(s)
    if s.built, return, end
    nb = numel(s.Bodies);
    % units: each rigid body, and each group of elastic bodies tied together
    par = 1:nb;
    for T = s.Ties
        a = T.master;  while par(a) ~= a, a = par(a); end
        b = T.slave;   while par(b) ~= b, b = par(b); end
        par(max(a, b)) = min(a, b);
    end
    top = zeros(1, nb);
    for i = 1:nb, a = i;  while par(a) ~= a, a = par(a); end, top(i) = a; end
    roots = unique(top, 'stable');  nu = numel(roots);
    s.bodyUnit = zeros(1, nb);
    s.bodyOffset = zeros(1, nb);
    off = 0;
    eoff = 0;
    s.anchors = struct('body', {}, 'frame', {}, 'x0', {}, 'Psi', {});
    s.units = [];
    for u = 1:nu
        members = find(top == roots(u));  s.bodyUnit(members) = u;
        if isscalar(members) && s.Bodies{members}.Rigid
            b = s.Bodies{members};  [mb, c, J] = b.massProperties();
            B = struct('elastic', false, 'm', mb, 'J', J, 'c0', c(:), 'Xl', b.X - c(:)', ...
                       'ne', 0, 'Phi', zeros(3*b.nNode, 0), 'Q', [], 's', [], 'CR', [], ...
                       'Kc', [], 'Dc', [], 'Mee', [], 'Psi', zeros(6, 0, 0), 'cpl', [], ...
                       'grounded', false, 'Q9', [], 'CRx', [], 'Qaz', [], 'Qju', [], 'Qjz', []);
        else
            noff = cumsum([0, cellfun(@(b) b.nNode, s.Bodies(members))]);
            s.bodyOffset(members) = noff(1:end-1);
            B = s.elasticGroup(members);
        end
        B.cols = off + (1:6 + B.ne);  B.ecols = eoff + (1:B.ne);
        off = off + 6 + B.ne;  eoff = eoff + B.ne;
        s.units = [s.units, B];
        s.anchors(u) = struct('body', u, 'frame', true, 'x0', zeros(3, 1), 'Psi', zeros(6, B.ne));
    end
    s.nv = off;
    q = struct('r', reshape([s.units.c0], 3, nu), 'R', repmat(eye(3), 1, 1, nu), 'e', zeros(eoff, 1));

    % reference points: rigid bodies -> unit anchor; couplings -> own
    % anchor (the first coupling of a free elastic group is its frame)
    s.refFrames = struct('anchor', {}, 's', {}, 'x0', {});
    for k = 1:numel(s.Refs)
        u = s.bodyUnit(s.Refs(k).body);
        x = s.Refs(k).x(:);
        B = s.units(u);
        if isempty(s.Refs(k).nodes)
            s.refFrames(k) = struct('anchor', u, 's', x - B.c0, 'x0', x);
        elseif ~B.grounded && B.cpl(1) == k
            s.refFrames(k) = struct('anchor', u, 's', zeros(3, 1), 'x0', x);
        else
            s.anchors(end+1) = struct('body', u, 'frame', false, 'x0', x - B.c0, ...
                                 'Psi', B.Psi(:,:,B.cpl == k));
            s.refFrames(k) = struct('anchor', numel(s.anchors), 's', zeros(3, 1), 'x0', x);
        end
    end

    % joints (zero at the defining configuration), supports of reference
    % points as joints to ground, supports of joints as locks
    s.jointFrames = struct('A', {}, 'B', {}, 'a', {}, 'b', {}, 'QA', {}, 'QB', {}, 'free', {}, 'k', {}, ...
                  'c', {}, 'pairs', {});
    for J = s.Joints
        s.addJoint(J.A, J.B, J.p(:), J.R, J.free, J.k, J.damp);
    end
    s.locks = struct('joint', {}, 'dir', {}, 'value', {});
    for f = s.Fixes
        dirs = f.dirs;
        if strcmp(f.kind, 'nodes'), continue, end % in the elastic group reduction
        if strcmp(f.kind, 'ref')
            if isempty(dirs), dirs = 1:6; end
            s.addJoint(0, f.idx, s.Refs(f.idx).x(:), eye(3), ~ismember(1:6, dirs), ...
                       zeros(6 - numel(unique(dirs)), 1), 0);
        else
            J = s.Joints(f.idx);  if isempty(dirs), dirs = 1:nnz(J.free); end
            for d = dirs(:)'
                v = 0;  if ~isempty(J.position), v = J.position(d); end
                s.locks(end+1) = struct('joint', f.idx, 'dir', d, 'value', v);
            end
        end
    end
    for u = find([s.units.grounded]) % elastic groups clamped on nodes
        s.jointFrames(end+1) = struct('A', 0, 'B', u, 'a', s.units(u).c0, 'b', zeros(3, 1), 'QA', eye(3), ...
                             'QB', eye(3), 'free', false(1, 6), 'k', zeros(6, 1), 'c', zeros(6, 1), ...
                             'pairs', rotationPairs(false(1, 6)));
    end
    s.transData = struct('terms', {}, 'coef', {}, 'k', {}, 'c', {});
    for T = s.Transmissions
        s.transData(end+1) = struct('terms', T.terms, 'coef', T.coef, 'k', T.k, 'c', T.damp);
    end

    % signals
    s.nodeSets = struct('body', {}, 'u0', {}, 'Phi', {}, 'x0', {});
    s.inSignals = s.compileSignals(s.Inputs);  s.outSignals = s.compileSignals(s.Outputs);

    % initial configuration: joint positions, transmissions unstretched,
    % elastic bodies undeformed; then the independent constraint rows
    s.keepRows = [];
    q = s.initialConfiguration(q);
    [~, G] = s.constraints(q, false);
    if ~isempty(G)
        [~, Rq, p] = qr(G', 'vector');
        r = sum(abs(diag(Rq)) > 1e-10*max(abs(diag(Rq))));
        s.keepRows = sort(p(1:r));
    end
    s.q0 = q;
    s.nu0 = s.initialVelocities(q);
    s.built = true;
end
