function B = reduceGroup(bodies, cpl, ties, sup, loads, nModes, name)
% Floating-frame data of an elastic group (bodies tied together). The frame
% is the first coupling (clamped), or, for a group supported on nodes,
% fixed at the first coupling or the centroid. Coupled nodes follow their
% coupling rigidly (RBE2), tie slaves the master surface; the other
% couplings are 6-dof interfaces. Craig-Bampton modes plus static shapes
% of the node-set loads; inertia invariants Q_ab = U_a' Ms U_b with
% U = [X0, Phi].
X = cell2mat(cellfun(@(b) b.X, bodies(:), 'UniformOutput', false));  nn = size(X, 1);
K = cell(1, numel(bodies));  M = K;  D = K;
for i = 1:numel(bodies)
    [K{i}, M{i}] = bodies{i}.matrices();
    D{i} = bodies{i}.Damping(1)*M{i} + bodies{i}.Damping(2)*K{i};
end
K = blkdiag(K{:});  M = blkdiag(M{:});  D = blkdiag(D{:});
for t = ties(arrayfun(@(t) ~isempty(t.k), ties))          % elastic interfaces: K += G' Ds G
    ns = numel(t.snodes);  kn = t.k(1);  kt = t.k(end);
    nnT = permute(t.n, [2 3 1]).*permute(t.n, [3 2 1]);
    Db = reshape(t.A, 1, 1, []).*(kt*eye(3) + (kn - kt)*nnT);
    [p, q, i] = ndgrid(1:3, 1:3, 1:ns);
    Ds = sparse(3*(i(:) - 1) + p(:), 3*(i(:) - 1) + q(:), Db(:));
    G = gapOperator(t, nn);
    K = K + G'*Ds*G;
end
grounded = ~isempty(sup);  nc = numel(cpl);
assert(grounded || nc > 0, 'fmbd:System', ['%s: an elastic body needs a coupling (s.couple), ' ...
       'which carries its frame, or node supports (s.fix).'], name);

% roles of the nodes: coupled, tie slave, supported dofs, free
owner = zeros(nn, 1);
for c = 1:nc
    if any(owner(cpl(c).nodes))
        error('fmbd:System', '%s: couplings %s and %s share nodes.', name, ...
              cpl(max(owner(cpl(c).nodes))).Name, cpl(c).Name);
    end
    owner(cpl(c).nodes) = c;
end
rigidTies = ties(arrayfun(@(t) isempty(t.k), ties));
slave = false(nn, 1);
for t = rigidTies, slave(t.snodes) = true; end
assert(~any(slave & owner > 0), 'fmbd:System', '%s: a tied node is also in a coupling.', name);
for t = rigidTies
    assert(~any(slave(t.mnodes(:))), 'fmbd:System', '%s: chained ties (a master node is tied).', name);
end
fixed = false(3*nn, 1);  fixed(sup) = true;
assert(~any(fixed(reshape(3*find(owner > 0 | slave)' - [2; 1; 0], [], 1))), 'fmbd:System', ...
       '%s: supported nodes cannot be coupled or tied.', name);
iface = 1:nc;  if ~grounded, iface = 2:nc; end          % the first coupling is the frame
freeDof = find(~fixed & ~reshape(repmat((owner > 0 | slave)', 3, 1), [], 1));
nf = numel(freeDof);  nq = nf + 6*numel(iface);
I = freeDof;  J = (1:nf)';  V = ones(nf, 1);
for ii = 1:numel(iface)                                 % u = u_c + theta_c x (x - x_c)
    c = iface(ii);  pn = cpl(c).nodes(:);  d = X(pn, :) - cpl(c).x;  base = nf + 6*(ii - 1);  np = numel(pn);
    r = 3*pn' - [2; 1; 0];
    I = [I; r(:); reshape(repmat(r, 3, 1), [], 1)];                                  %#ok<AGROW>
    J = [J; repmat(base + (1:3)', np, 1); reshape(repmat(base + [4 5 6], 3, np), [], 1)]; %#ok<AGROW>
    H = zeros(3, 3, np);                                % -hat(d) per node
    H(1,2,:) = d(:,3);  H(1,3,:) = -d(:,2);  H(2,1,:) = -d(:,3);
    H(2,3,:) = d(:,1);  H(3,1,:) = d(:,2);   H(3,2,:) = -d(:,1);
    V = [V; ones(3*np, 1); H(:)];                                                      %#ok<AGROW>
end
T = sparse(I, J, V, 3*nn, nq);
if ~isempty(rigidTies)                                  % tie slaves follow the master surface
    P = sparse(3*nn, 3*nn);
    for t = rigidTies, P = P - gapOperator(t, nn, true); end
    T = T + P*T;
end
Kr = T'*K*T;  Mr = T'*M*T;  Dr = T'*D*T;  Kr = (Kr + Kr')/2;  Mr = (Mr + Mr')/2;
ib = (nf + 1:nq)';
[~, ~, Vr] = fmbd.craigBampton(Kr, Mr, ib, nModes);
if ~isempty(loads)                                      % residual attachment shapes of the loads
    in = (1:nf)';  fr = T'*loads;  xs = zeros(nq, size(loads, 2));
    xs(in, :) = decomposition(Kr(in, in), 'chol')\full(fr(in, :));
    xs = xs - Vr*((Vr'*Mr*Vr)\(Vr'*Mr*xs));
    [Qs, Ss] = eig(full(xs'*Mr*xs + (xs'*Mr*xs)')/2, 'vector');
    kp = Ss > 1e-10*max([Ss; eps]);
    Vr = [Vr, xs*Qs(:, kp)./sqrt(Ss(kp))'];
end
ne = size(Vr, 2);
Kc = Vr'*Kr*Vr;  Mc = Vr'*Mr*Vr;  Dc = Vr'*Dr*Vr;
Phi = full(T*Vr);
Psi = zeros(6, ne, nc);
for ii = 1:numel(iface), Psi(:,:,iface(ii)) = Vr(nf + 6*(ii - 1) + (1:6), :); end
if nc > 0, x0 = cpl(1).x(:); else, x0 = mean(X, 1)'; end
Xl = X - x0';
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
B = struct('elastic', true, 'm', full(sum(Ms(:))), 'J', [], 'c0', x0, 'Xl', Xl, 'ne', ne, ...
           'Phi', Phi, 'Q', Q, 's', s, 'CR', CR, 'Kc', (Kc + Kc')/2, 'Dc', (Dc + Dc')/2, ...
           'Mee', (Mc + Mc')/2, 'Psi', Psi, 'cpl', [], 'grounded', grounded);
B = addReshapes(B);
end
