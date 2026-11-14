function x = perturb(s, x, dq)
%PERTURB Move a state along a tangent direction dq.
    build(s);
    [q, nu] = startState(s, x);  u = [];
    if isfield(x, 'u') && ~isempty(x.u), u = x.u(:, end); end
    x = state(s, retract(s, q, dq), nu, inputVector(s, u));
end
