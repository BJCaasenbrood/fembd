function u = inputVector(s, u)
% Input signals as a column (default zero).
    if isempty(u), u = zeros(numel(s.inSignals), 1); end
    u = u(:);
    assert(numel(u) == numel(s.inSignals), 'fmbd:System', '%d input value(s) expected.', numel(s.inSignals));
end
