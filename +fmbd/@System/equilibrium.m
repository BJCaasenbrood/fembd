function [x, info] = equilibrium(s, opts)
%EQUILIBRIUM Static equilibrium under gravity, springs and constant
%   inputs: x = s.equilibrium(Input = u0) starts from the initial
%   configuration (or Initial) and returns the state (as from simulate,
%   one step, at rest, with the inputs in x.u). Coordinates without
%   stiffness (free mechanisms) stay where they are. Use it with
%   simulate(Initial = x), linearize, statespace, modes and plot.
%
%   Newton iterations with a backtracking line search. Converged when
%   the force residual is below Tol (relative to the forces) and the
%   constraint residual below Tol (relative to the model size). If the
%   residual stops decreasing, it is accepted at round-off level, and
%   otherwise the iterations stop with a warning (e.g. a free mechanism
%   with a net force has no equilibrium).
%   [x, info] = ... also returns info: converged, iterations, force
%   (relative force residual) and constraint (constraint residual).
%
%   Options: Input, Initial, Tol (default 1e-10), MaxIter (default 100)
    arguments
        s
        opts.Input double = []
        opts.Initial struct = struct([])
        opts.Tol (1,1) double {mustBePositive} = 1e-10
        opts.MaxIter (1,1) double = 100
    end
    s.build();
    [q, ~] = s.startState(opts.Initial);
    n = s.nv;
    nu = zeros(n, 1);
    u = s.inputVector(opts.Input);
    L = max([1; abs(q.r(:))]); % model size, for the constraint tolerance
    F = s.forces(q, nu, u, false);  [g, G] = s.constraints(q);
    lam = lsqminnorm(G', F);
    r = [F - G'*lam; g]; % force balance, constraints
    err = @(r, F) [norm(r(1:n))/max(norm(F), 1), max(abs([r(n+1:end); 0]))/L];
    merit = @(e) e(1) + e(2);
    best = inf;
    stall = 0;
    conv = false;
    Kt = [];
    age = 0;
    for it = 1:opts.MaxIter
        e = err(r, F);
        if all(e < opts.Tol), conv = true;  break, end
        m = max(e/opts.Tol); % progress measure
        if m < 0.5*best
            best = m;  stall = 0;
        elseif age > 1
            Kt = []; % slow with an old tangent: refresh it first
        else
            stall = stall + 1; % slow even with a fresh tangent
        end
        if stall >= 3 % plateau: round-off floor, or no equilibrium
            conv = all(e < 1e3*opts.Tol);
            break
        end
        if isempty(Kt) || age >= 10 % modified Newton: reuse the tangent
            Kt = s.frameStiffness(q, u, lam);  age = 0;
        end
        % Newton step split by the constraints: correct them exactly, then
        % balance the forces in the directions they allow (null space of G)
        r1 = r(1:n);
        dqc = -lsqminnorm(G, g);
        Nn = null(G);  if isempty(G), Nn = eye(n); end
        z  = lsqminnorm(Nn'*Kt*Nn, Nn'*(r1 - Kt*dqc)); % free mechanisms stay
        dq = dqc + Nn*z;
        dlam = lsqminnorm(G', r1 - Kt*dq);
        a = 1;  m0 = merit(e);
        for ls = 1:6 % backtracking line search
            q1 = s.retract(q, a*dq);  lam1 = lam + a*dlam;
            F1 = s.forces(q1, nu, u, false);  [g1, G1] = s.constraints(q1);
            r1n = [F1 - G1'*lam1; g1];
            if merit(err(r1n, F1)) < m0 || ls == 6, break, end
            a = a/2;
        end
        [q, lam, F, G, g, r] = deal(q1, lam1, F1, G1, g1, r1n);
        age = age + 1;
    end
    e = err(r, F);
    info = struct('converged', conv, 'iterations', it, 'force', e(1), 'constraint', e(2)*L);
    if ~conv
        warning('fmbd:equilibrium', ['No equilibrium found after %d iterations (force residual ' ...
            '%.1e, constraint residual %.1e); a free coordinate may carry a net force.'], ...
            it, e(1), e(2)*L);
    end
    x = s.state(q, nu, u);
end
