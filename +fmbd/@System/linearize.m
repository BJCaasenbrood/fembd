function lin = linearize(s, x, opts)
%LINEARIZE Linear model about a state at rest (default: the initial
%   configuration), e.g. from equilibrium:
%     M z'' + C z' + K z = B u,   y = Cz z + Cv z' + Ca z''
%   in independent coordinates z (dq = T*z satisfies the linearised
%   constraints, T orthonormal). K includes the stiffness from gravity,
%   constraint forces, preloaded springs and the operating-point inputs
%   (they act through the geometry, e.g. a torque holding a spring).
%   Returns a struct with M, C, K, B, Cz, Cv, Ca, T, InputName and
%   OutputName. Option Input = u0 (default: x.u, e.g. from equilibrium).
    arguments
        s
        x struct = struct([])
        opts.Input double = []
    end
    s.build();
    [q, nu] = s.startState(x);
    if any(nu)
        warning('fmbd:linearize', 'The state is moving; the model is linearised at its velocities.');
    end
    n = s.N;  u0 = opts.Input;              % the inputs at the operating point
    if isempty(u0) && ~isempty(x) && isfield(x, 'u') && ~isempty(x.u), u0 = x.u(:, end); end
    u0 = s.inputVector(u0);
    [~, G] = s.constraints(q);
    T = null(G);  if isempty(G), T = eye(n); end
    [F, ~, ~, Mq] = s.forces(q, nu, u0);  nc = size(G, 1);
    sol = lsqminnorm([Mq, G'; G, zeros(nc)], [F; -s.convective(q, nu)]);   % constraint forces
    lam = sol(n+1:end);                         % (dynamic, so also away from equilibrium)
    Kq = s.stiffness(q, u0, lam, nu);
    Cq = zeros(n);  e = 1;                      % forces are at most quadratic in nu,
    for j = 1:n                                 % so central differences are exact
        d = zeros(n, 1);  d(j) = e;
        Cq(:, j) = -(s.forces(q, nu + d, u0) - s.forces(q, nu - d, u0))/(2*e);
    end
    Bq = zeros(n, numel(s.In));
    for i = 1:numel(s.In)                       % forces are linear in the inputs
        du = zeros(numel(s.In), 1);  du(i) = 1;
        Bq(:, i) = s.forces(q, nu, u0 + du) - s.forces(q, nu, u0);
    end
    no = numel(s.Out);  Yq = zeros(no, n);  Yv = zeros(no, n);  Ya = zeros(no, n);
    for i = 1:no
        sg = s.Out(i);  [~, Jc] = s.coord(sg.kind, sg.idx, q);
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
