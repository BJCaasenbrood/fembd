function lin = linearize(s, x, opts)
%LINEARIZE Linear model about a state at rest.
%   lin = s.linearize(x) returns M, C, K, B, Cz, Cv, Ca and T of
%     M z'' + C z' + K z = B u,   y = Cz z + Cv z' + Ca z''
%   in independent coordinates z, with dq = T*z. K includes gravity,
%   constraint forces and preload. Default x: the initial configuration.
%   Option Input (default: x.u).
    arguments
        s
        x struct = struct([])
        opts.Input double = []
    end
    build(s);
    [q, nu] = startState(s, x);
    if any(nu)
        warning('fmbd:linearize', 'The state is moving; the model is linearised at its velocities.');
    end
    n = s.nv;  u0 = opts.Input; % the inputs at the operating point
    if isempty(u0) && ~isempty(x) && isfield(x, 'u') && ~isempty(x.u), u0 = x.u(:, end); end
    u0 = inputVector(s, u0);
    [~, G] = constraints(s, q);
    T = null(G);  if isempty(G), T = eye(n); end
    [F, ~, ~, Mq] = forces(s, q, nu, u0);  nc = size(G, 1);
    sol = lsqminnorm([Mq, G'; G, zeros(nc)], [F; -convective(s, q, nu)]); % constraint forces
    lam = sol(n+1:end); % (dynamic, so also away from equilibrium)
    Kq = stiffness(s, q, u0, lam, nu);
    Cq = zeros(n);  e = 1; % exact: the forces are quadratic in nu
    for j = 1:n
        d = zeros(n, 1);  d(j) = e;
        Cq(:, j) = -(forces(s, q, nu + d, u0) - forces(s, q, nu - d, u0))/(2*e);
    end
    Bq = zeros(n, numel(s.inSignals));
    for i = 1:numel(s.inSignals) % forces are linear in the inputs
        du = zeros(numel(s.inSignals), 1);  du(i) = 1;
        Bq(:, i) = forces(s, q, nu, u0 + du) - forces(s, q, nu, u0);
    end
    no = numel(s.outSignals);
    Yq = zeros(no, n);
    Yv = zeros(no, n);
    Ya = zeros(no, n);
    for i = 1:no
        sg = s.outSignals(i);  [~, Jc] = coord(s, sg.kind, sg.idx, q);
        switch sg.type
            case 'displacement', Yq(i,:) = Jc(sg.dir, :);
            case 'velocity',     Yv(i,:) = Jc(sg.dir, :);
            otherwise,           Ya(i,:) = Jc(sg.dir, :);
        end
    end
    lin = struct('M', T'*Mq*T, 'C', T'*Cq*T, 'K', T'*Kq*T, 'B', T'*Bq, 'Cz', Yq*T, ...
                 'Cv', Yv*T, 'Ca', Ya*T, 'T', T, 'InputName', {s.InputName}, ...
                 'OutputName', {s.OutputName});
    lin.M = (lin.M + lin.M')/2;
end
