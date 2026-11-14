function k = refIndex(s, id)
% Reference point index of id; 0 for ground.
    k = 0;
    if (ischar(id) || isstring(id)) && strcmpi(id, 'ground'), return, end
    [kind, k] = lookup(s, id, true);
    assert(strcmp(kind, 'ref'), 'fmbd:joint', ...
        'Joint ends must be rigid bodies, couplings (s.couple) or ''ground''.');
end
