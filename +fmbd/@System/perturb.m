function x = perturb(s, x, dq)
%PERTURB Move a state along a tangent direction dq.
    s.build();
    [q, nu] = s.startState(x);  u = [];
    if isfield(x, 'u') && ~isempty(x.u), u = x.u(:, end); end
    x = s.state(s.retract(q, dq), nu, s.inputVector(u));
end
