function u = inputVector(s, u)
% Input signals as a column (default zero).
    if isempty(u), u = zeros(numel(s.In), 1); end
    u = u(:);
    assert(numel(u) == numel(s.In), 'fmbd:System', '%d input value(s) expected.', numel(s.In));
end
