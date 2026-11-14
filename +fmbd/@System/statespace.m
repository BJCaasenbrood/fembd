function sys = statespace(s, x, opts)
%STATESPACE Modal state-space model.
%   sys = s.statespace(x) linearises about x and returns A, B, C, D in
%   modal coordinates [eta; deta], with signal names and frequencies f.
%   Options:
%     Modes             number of modes (default: all)
%     Damping           modal damping ratio, added to the physical dampers
%     StaticCorrection  keep the static gain when modes are cut (default true)
%     Input             inputs at the operating point (default: x.u)
%     Format            "struct" (default) or "ss" (Control System Toolbox)
    arguments
        s
        x struct = struct([])
        opts.Modes (1,1) double {mustBePositive} = inf
        opts.Damping double = 0
        opts.StaticCorrection (1,1) logical = true
        opts.Input double = []
        opts.Format (1,1) string {mustBeMember(opts.Format, ["struct", "ss"])} = "struct"
    end
    lin = s.linearize(x, Input = opts.Input);
    K = (lin.K + lin.K')/2;
    M = lin.M;
    B = lin.B;
    nz = size(K, 1);
    [P, lam] = eig(K, M, 'vector');
    [lam, i] = sort(real(lam));  P = real(P(:, i));
    lam(abs(lam) < 1e-10*max(abs(lam))) = 0; % round-off rigid/mechanism modes: exactly 0
    nm = min(opts.Modes, nz);  V = P(:, 1:nm);
    if opts.StaticCorrection && nm < nz && ~isempty(B)
        sigma = 0;  if any(lam(1:nm) <= 0), sigma = -1e-6*max(abs(lam)); end
        X = (K - sigma*M)\B;  X = X - V*(V'*(M*X));
        [Q, S] = eig((X'*M*X + (X'*M*X)')/2, 'vector');  kp = S > 1e-12*max([S; eps]);
        W = [V, X*Q(:,kp)./sqrt(S(kp))'];
        [P2, l2] = eig((W'*K*W + (W'*K*W)')/2, (W'*M*W + (W'*M*W)')/2, 'vector');
        [l2, i] = sort(real(l2));
        V = W*real(P2(:, i));
        lam = l2;
        lam(abs(lam) < 1e-10*max(abs(lam))) = 0;
    else
        lam = lam(1:nm);
    end
    V = V./sqrt(dot(V, M*V));  nr = size(V, 2);
    zeta = opts.Damping(:);  if isscalar(zeta), zeta = zeta*ones(nr, 1); end
    w = sqrt(max(lam, 0));
    Z = diag(2*zeta.*w) + V'*lin.C*V;  Om2 = V'*lin.K*V;
    A = [zeros(nr), eye(nr); -Om2, -Z];
    Bs = [zeros(nr, size(B, 2)); V'*B];
    Ca = lin.Ca*V;  Bm = V'*B; % accelerations: eta'' = -Om2 eta - Z eta' + Bm u
    Cs = [lin.Cz*V - Ca*Om2, lin.Cv*V - Ca*Z];  Ds = Ca*Bm;
    sys = struct('A', A, 'B', Bs, 'C', Cs, 'D', Ds, 'InputName', {lin.InputName}, ...
                 'OutputName', {lin.OutputName}, 'f', sign(lam).*sqrt(abs(lam))/(2*pi));
    if opts.Format == "ss"
        sys = ss(A, Bs, Cs, Ds, 'InputName', lin.InputName, 'OutputName', lin.OutputName);
    end
end
