function x = perturb(s, x, dq)
%PERTURB State x moved along the tangent direction dq (nDof x 1), e.g.
%   a mode shape from modes: s.plot(s.perturb(x, 0.05*V(:,k))).
    s.build();
    [q, nu] = s.startState(x);  u = [];
    if isfield(x, 'u') && ~isempty(x.u), u = x.u(:, end); end
    x = s.state(s.retract(q, dq), nu, s.inputVector(u));
end
