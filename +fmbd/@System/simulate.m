function [t, y, x] = simulate(s, tspan, opts)
%SIMULATE Simulate the system in time.
%   [t, y, x] = s.simulate([t0 t1]) starts in the initial configuration and
%   returns the time t, the outputs y (struct) and the states x, with
%   fields r, R (frame positions and rotations per unit), e (modal
%   coordinates), nu (velocities) and u (inputs).
%   Options:
%     Step     time step (default 1e-3)
%     Input    controller u = f(t, y), using the outputs of the previous step
%     Initial  start state, e.g. from equilibrium
%     Rho      numerical damping, 0..1 (default 0.8; 1 = none)
%     Tol      Newton tolerance (default 1e-8)

   arguments
        s
        tspan (1,2) double
        opts.Step (1,1) double {mustBePositive} = 1e-3
        opts.Input = []
        opts.Rho (1,1) double {mustBeInRange(opts.Rho, 0, 1)} = 0.8
        opts.Tol (1,1) double {mustBePositive} = 1e-8
        opts.Initial struct = struct([])
    end
    build(s);
    nt = max(1, round(diff(tspan)/opts.Step));  h = diff(tspan)/nt;
    t = tspan(1) + h*(0:nt)';
    r = opts.Rho;
    am = (2*r - 1)/(r + 1);
    af = r/(r + 1);
    gm = 0.5 + af - am;
    be = 0.25*(gm + 0.5)^2;
    ca = (1 - af)/(1 - am);
    nb = numel(s.units);  n = s.nv;

    [q, nu] = startState(s, opts.Initial);  u = zeros(numel(s.inSignals), 1);
    [F, ~, ~, Mq] = forces(s, q, nu, u);
    [~, G] = constraints(s, q);  nc = size(G, 1);
    sol = [Mq, G'; G, zeros(nc)]\[F; -convective(s, q, nu)];
    wd = sol(1:n);
    lam = sol(n+1:end);
    acc = wd;
    y0 = outputs(s, q, nu, wd);
    if ~isempty(opts.Input) % first input, then consistent accelerations
        u = opts.Input(t(1), y0);  F = forces(s, q, nu, u);
        sol = [Mq, G'; G, zeros(nc)]\[F; -convective(s, q, nu)];
        wd = sol(1:n);  lam = sol(n+1:end);  acc = wd;  y0 = outputs(s, q, nu, wd);
    end

    x = struct('r', zeros(3, nb, nt+1), 'R', zeros(3, 3, nb, nt+1), ...
               'e', zeros(numel(q.e), nt+1), 'nu', zeros(n, nt+1), 'u', zeros(numel(s.inSignals), nt+1));
    Y = zeros(numel(s.outSignals), nt+1);
    x.r(:,:,1) = q.r;
    x.R(:,:,:,1) = q.R;
    x.e(:,1) = q.e;
    x.nu(:,1) = nu;
    x.u(:,1) = u;
    Y(:,1) = cell2mat(struct2cell(y0));
    yk = y0;
    for k = 1:nt
        if ~isempty(opts.Input), u = opts.Input(t(k), yk); end
        w = wd;  prev = inf;
        for it = 1:30
            an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
            q1 = retract(s, q, h*(nu + h*(0.5 - be)*acc + h*be*an));
            nu1 = nu + h*(1 - gm)*acc + h*gm*an;
            if it == 1 % tangents once per step (modified Newton)
                [F, Kt, Ct, Mq] = forces(s, q1, nu1, u);
            else
                [F, ~, ~, Mq] = forces(s, q1, nu1, u, false);
            end
            [g, G] = constraints(s, q1);
            res = [Mq*w - F + G'*lam; g/(h^2*be*ca)];
            S = [Mq + h*gm*ca*Ct + h^2*be*ca*Kt, G'; G, zeros(nc)];
            dx = -S\res;
            w = w + dx(1:n);  lam = lam + dx(n+1:end);
            e = norm(dx(1:n))/(1 + norm(w));
            if e <= opts.Tol || (e >= prev && e < 1e3*opts.Tol), break, end % converged, or at round-off
            prev = e;
        end
        if it == 30 && e > 1e3*opts.Tol
            warning('fmbd:simulate', 'Newton did not converge at t = %g s; try a smaller Step.', t(k+1));
        end
        an = ((1 - af)*w + af*wd - am*acc)/(1 - am);
        q = retract(s, q, h*(nu + h*(0.5 - be)*acc + h*be*an));
        nu = nu + h*(1 - gm)*acc + h*gm*an;
        if nc > 0 % project onto G nu = 0 (mass-orthogonal)
            [~, G] = constraints(s, q);  Mq = massMatrix(s, q);
            p = [Mq, G'; G, zeros(nc)]\[zeros(n, 1); -G*nu];  nu = nu + p(1:n);
        end
        acc = an;  wd = w;
        yk = outputs(s, q, nu, w);
        x.r(:,:,k+1) = q.r;
        x.R(:,:,:,k+1) = q.R;
        x.e(:,k+1) = q.e;
        x.nu(:,k+1) = nu;
        x.u(:,k+1) = u;
        Y(:,k+1) = cell2mat(struct2cell(yk));
    end
    y = struct();
    for i = 1:numel(s.outSignals)
        y.(matlab.lang.makeValidName(s.Outputs(i).Name)) = Y(i,:)';
    end
end
